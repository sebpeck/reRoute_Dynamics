'''
Geography_Tools.py
'''
import pandas as pd
import geopandas as gpd
import shapely
import numpy as np
import datetime as dt
from alive_progress import alive_bar
import inspect
import rasterio
import os
import scipy
import pyogrio as pio
import json
import RouteMap as rm
import random

def haversine_formula(x1, y1, x2, y2):
    '''
    haversine_formula takex the latitude and longitude of two separate points,
    and calculates the distance in kilometers between those two points as a crow flies.
    
    Params:
    x1 - latitude pt 1 in degrees
    y1 - longitude pt 1 in degrees
    x2 - latitude pt 2 in degrees
    y2 - longitude pt 2 in degrees
    
    Returns:
    distance between the points in kilometers.
    
    Notes:
    11/20/2024 - Efficiency can be improved by having radian conversion be done externally.
    '''
    
    # convert the degree coordinates to radians
    lat1=np.radians(x1)
    lon1=np.radians(y1)
    lat2=np.radians(x2)
    lon2=np.radians(y2)
    
    # Get the radius of earth in kilometers
    R = 6373.0
    
    # calculate the difference in radians
    dlon=lon2-lon1
    dlat=lat2-lat1
    
    # calculate the distance using the haversine formula.
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distance = R * c
    
    
    # return the distance in km.
    if distance < .0001:
        distance=.0001
    return distance


def compass_heading(bearing):
    '''
    compass_heading converts a bearing in degrees to a
    compass value like North, South, or East.
    
    Params:
    bearing - the bearing value, in degrees, as float.
    
    Returns:
    string representation of compass heading.
    
    Notes: 
    11/18/2024 - There is a more mathematically elegant solution to this, but I'm not bothering
                 with it right now.
    '''
    
    # set up a list of possible compass directions
    possible_dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"]
    
    # get the index of the 45 quadrants the angle bearing is in
    bearing_index = int(np.round(bearing/45))
    
    # make sure the index isn't negative
    if bearing_index < 0:
        bearing_index=bearing_index+8
        
    # return the corresponding compass direction
    return possible_dirs[bearing_index]


def heading_to_angle(heading):
    '''
    heading_to_angle() takes a compass heading of 8 directions,
    and converts it to an angle in degrees. 
    
    Params:
    heading - a string representing heading.
    
    Returns: a compass bearing angle as float.
    '''
    options = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    idx = options.index(heading)
    return (idx/len(options)*360)
    
    

def point_bearing(x1, y1, x2, y2, bearing_type = "Angle"):
    '''
    point_bearing takex in the latitude and longitude of two separate ponts,
    and calculates the bearing when travelling from point 1 and point 2.
    
    Params:
    x1 - latitude pt 1 in degrees
    y1 - longitude pt 1 in degrees
    x2 - latitude pt 2 in degrees
    y2 - longitude pt 2 in degrees
    bearing_type - determine the bearing output.
                    Options: Angle - return the angle in degrees
                             Compass - return the angle as a directional bearing (eg: N, E, S, W)
    
    Returns:
    bearing of the vector between the two points. 
    
    '''
    # convert to radians
    lat1=np.radians(x1)
    lon1=np.radians(y1)
    lat2=np.radians(x2)
    lon2=np.radians(y2)
    
    # run the bearing calculation
    x = np.cos(lat2) * np.sin(lon2-lon1)
    y = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2-lon1)
    bearing = np.arctan2(x, y)
    bearing = np.degrees(bearing)
    
    # check the return type
    if (bearing_type == "Angle"):
        return bearing
    elif (bearing_type == "Compass"):
        return compass_heading(bearing)


def get_bounding_box(shape):
    '''
    get_bounding_box() takes in a shapely shape,
    and then generates a polygon that is a rectangular bounding box.
    
    Parameters:
    shape: any shapely shape
    
    Returns: 
    shapely.polygon that fits around the provided shape
    '''
    # get the bound coordinates
    bbox_cords = shape.bounds
    
    # create a rectangle of those coordinates
    boundary_coords = ((bbox_cords[0], bbox_cords[1]),
                 (bbox_cords[0], bbox_cords[3]),
                 (bbox_cords[2], bbox_cords[3]),
                 (bbox_cords[2], bbox_cords[1]))
    
    # return a shapely polygon of those coordinates
    return shapely.Polygon(boundary_coords)


def interpolate_points(point_1, point_2, max_dist=1):
    '''
    interpolate_points takes two shapely points and a maximum distance,
    and then interpolates between the points if the maximum distance is exceeded.
    
    Params:
    point_1 - starting point as a shapely point of lat, lon
    point_2 - final point as a shapely point of lat, lon
    max_dist - maximum distance the points can be without interpolation. Default int of 1.
    
    Returns:
    an iterable of lon,lat shapely points. 
    '''
    
    # get distance in meters
    distance = haversine_formula(point_1.y, point_1.x, point_2.y, point_2.x)*1000
    if distance < max_dist:
        return [point_1]
    else:
        # calculate the number of points between to be interpolated
        num_interp = int(np.ceil(distance/max_dist)) - 1
        
        # calculate the distance between each interpolated point
        dx = distance/num_interp
        
        # generate a linestring
        line = shapely.LineString([point_1, point_2])
        
        # create the points list with the initial point
        points = [point_1]
        
        # loop through each distance from the origin 
        for dis in np.arange(dx/distance, 1, dx/distance):
            
            # interpolate a point
            i_point = line.interpolate(dis, normalized=True).coords[:][0]
            
            # convert it to a shapely point
            i_point = shapely.Point(i_point[0], i_point[1])
            
            # append the point
            points.append(i_point)
            
        # return the points
        return points
    
    
def repeat_id_remover(sequence):
    '''
    repeat_id_remover takes an iterable of IDs where -1 is invalid,
    and swaps out any repeat index with an invalid.
    
    Params:
    sequence - a sequence of id values in an iterable.
    
    Returns:
    sequence with repeated values swapped for a -1.
    '''
    
    # start the sequence with invalid
    sequence_value = -1
    
    # build a list for the new sequence
    new_sequence = []
    
    # loop through the old
    for item in sequence:
        
        # check if the item is the same as the saved val
        if (item == sequence_value):
            
            # if so, append invalid
            new_sequence.append(-1)
        else:
            
            # otherwise, append the item and swap the saved val for the item
            new_sequence.append(item)
            sequence_value=item
            
    # return the new sequence
    return new_sequence


### ---- SERIES QUERY METHODS ---- ###

def verbose_line_updater(message, reset=False):
    '''
    verbose_line_updater generates a verbose message with timestamps and
    method calls based on a passed message. Can reset printline if specified.
    
    Params:
    message - message, as str, to be printed with the line.
    reset - boolean to reset carrige or not. Default false.
    
    Returns:
    The string to be printed.
    '''
    
    # generate the string.
    data_string = "[{} -- {}.{}()] {}                  ".format(dt.datetime.now().time(),
                                                                "Geography_Tools",
                                                                str(inspect.stack()[1][3]),
                                                                message)
    # check if the carrage should be reset.
    if reset:
        print(data_string, end='\r')
    else:
        print(data_string)
        
    # return the carrage.
    return data_string


def query_elevation_changes(elev):
    '''
    query_elevation_change takes an iterable of elevation,
    and returns the elevation change between the point at the
    current index, and the next poinnt. 
    
    Params:
    elev - iterable of elevation data
    
    Returns:
    list of elevation changes to reach the next point.
    '''
    
    # convert to series
    elev_ser = pd.Series(elev)
    
    # shift the series by one
    shifted = elev_ser.shift(-1)
    
    # get the original final point
    final_point = elev[-1]
    
    # set the final point to replace the nan
    shifted.iloc[-1] = final_point
    
    # calculate the shift between the elevations and return.
    return list(elev_ser-shifted)
    
    
def query_stops(geometry, stop_table_path, key='stop_id', epsg_from = 4326, epsg_to=4326, margin=10, verbose=False):
    '''
    query_stops takes in a geometry point series for a given route, and a path to
    all stop ids and geometry, and then returns a list of stop id's in the order the
    bus will arrive.
    
    Params:
    geometry - a series of shapely points representing longitude and latitude
    stop_table_path - a path to the dataset containing stop geodata csv, as str
    key - column identifier for the stop id (or whatever value you want) as str
    epsg_from - the epsg projection the csv is in as an int, default is 4326
    epsg_to - the epsg projection the data should be reporojected to as an int, default is 4326
    margin - distance, in meters, that the stop is allowed to be from a given point to qualify, as int.
    verbose - boolean to enable verbosity. default False.
    
    Returns:
    list of stop ids (or other queried values from the data), with invalid points being labeled -1.
    '''
    
    # Check the shape boundaries for the geometric series
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))
    
    # read the data.
    if verbose: verbose_line_updater("Loading CSV: {}".format(stop_table_path))
    stop_table = pd.read_csv(stop_table_path)
    if verbose: verbose_line_updater("Loaded CSV.")
    
    # convert the data to a shapely geometry pandas series
    if verbose: verbose_line_updater("Loading stop geometry...")
    stop_table['geometry'] = stop_table.apply(lambda x: shapely.Point(x.stop_lon, x.stop_lat), axis=1)
    if verbose: verbose_line_updater("Loaded stop geometry.")
    
    # convert the data to a geodataframe, and re-set the crs to the specified
    if verbose: verbose_line_updater("Setting geodata: epsg{}-->epsg:{}".format(epsg_from, epsg_to))
    stop_table = gpd.GeoDataFrame(stop_table).set_crs(epsg = epsg_from).to_crs(epsg = epsg_to)
    if verbose: verbose_line_updater("Geodata Set.")
    
    # filter the stop table to only be the stops within the boundary
    if verbose: verbose_line_updater("Filtering stops...")
    stop_table = stop_table[stop_table['geometry'].apply(lambda x: shape_bounds.contains(x)) == True]
    if verbose: verbose_line_updater("Stops filtered.")
    
    
    # check for stops within a set distance of each geometry point
    stop_id_list = []
    
    if verbose: verbose_line_updater("Processing stops...")
    # loop through the geimetry
    
    with alive_bar(len(geometry)) as bar:
        for point in list(geometry):

            # re-set the pont name because i dont want to re-change them all
            route_pt = point

            # get the distances using the haversine formula, converted to meters
            stop_table['dist'] = stop_table.geometry.apply(lambda x: haversine_formula(x.y, x.x, route_pt.y, route_pt.x)*1000)

            # filter the data to yeild the closest stop
            closest_stop = stop_table[stop_table['dist'] == stop_table['dist'].min()].reset_index().iloc[0]

            # if the closest stop is within the margin, append it to the list. Otherwise, append -1.
            stop_id_list.append(int(closest_stop['dist']<margin)*closest_stop[key] - int(closest_stop['dist']>margin))
            bar()
    if verbose: verbose_line_updater("Stops processed. Returning values.")
    
        
    # return the stop id list with repeats removed, and the last index guaranteed to be a stop.
    return repeat_id_remover(stop_id_list)


def query_distance_traveled(geometry_series, verbose=False):
    '''
    query_distance_traveled takes in a series of points,
    and then calculates the distance traveled from point to point, in kilometers.

    Parameters:
    geometry_series - a series of shapely points of lat and lon.
    verbose - boolean to enable verbosity. default False.

    Returns:
    an iterable representing the change in distance between each point and the next.
    '''

    # get the current geometry series
    if verbose: verbose_line_updater("Loading geometry.")
    currents = geometry_series
    
    # shift the geometry series up by one, which will be the next ponts.
    if verbose: verbose_line_updater("Shifting geometry.")
    nexts = geometry_series.shift(-1)

    # combine the two into a dataframe.
    if verbose: verbose_line_updater("Combining geometry.")
    data = pd.concat([currents, nexts], axis=1, keys=['current', 'next']).reset_index(drop=True)

    # get the final point and make sure it's the last value.
    final_point = list(data['current'])[-1]

    # set the nan value in the beginning of the 'nexts' to the initial, representing no distance traveled.
    data.iloc[-1, data.columns.get_loc("next")] = final_point

    # apply the haversine formula to the data, which will yeild a series of changes in distance.
    if verbose: verbose_line_updater("Calculating distances.")
    dXs = data.apply(lambda x: haversine_formula(x.current.y, x.current.x, x.next.y, x.next.x), axis=1)
    
    if verbose: verbose_line_updater("Point distances processed. Returning values.")

    # return the series of the changes in distance, in km, make sure it never doesn't travel.
    return dXs.values


def query_signals(geometry, signal_data_path, key="SIGNAL_ID", epsg = 4326, margin = 10, verbose=False):
    '''
    query_signals takes in a geometric series of shapely points,
    a path to stoplight signal data, and returns a list of the signal id at
    each point in that series.
    
    Parameters:
    geometry - a series of shapely points of lat and lon.
    signal_data_path - path, as str, to a shapefile containing stop signal data.
    key - the column identifier of signal id (or other target value) to be queried
    epsg - the epsg the data will be re-set as. Default is 4326.
    margin - the distance, in meters, that the signal is allowed to be from a given point to qualify.
    verbose - boolean to enable verbosity. Default False.
    
    Returns:
    list of signal IDs, or other info about the signal, as specified.
    -1 if none found.
    '''
    
    # get the bounding box of the geometry
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))
    
    # Get the signal data.
    if verbose: verbose_line_updater("Loading shapefile & setting EPSG: {}".format(signal_data_path))
    signals = pio.read_dataframe(signal_data_path).to_crs(epsg = epsg)
    if verbose: verbose_line_updater("Shapefile laded with EPSG:{}.".format(epsg))
    
    # get the signals that are within the bounds of the data.
    if verbose: verbose_line_updater("Filtering sinals...")
    signals = signals[signals.apply(lambda x: shape_bounds.contains(x.geometry), axis=1).reset_index(drop=True) == True]
    if verbose: verbose_line_updater("Signals filtered.")
    
    # Generate empty list for signal ids
    signal_id_list = []
    
    # loop through the geometry
    if verbose: verbose_line_updater("Processing signals...")
    with alive_bar(len(geometry)) as bar:
        for point in list(geometry):

            # use the haversine formula to get the distances between the current point and each point in the series, in meters
            signals['dist'] = signals.geometry.apply(lambda x: haversine_formula(x.y, x.x, point.y, point.x)*1000)

            # get the closest signal through the minimum distance
            closest_signal = signals[signals['dist'] == signals['dist'].min()].reset_index().iloc[0]

            # check if the signal is within the margin, and if so, append it to the list. If not, append -1.
            signal_id_list.append((int(closest_signal['dist'])<margin)*closest_signal[key] - int(closest_signal['dist']>margin))
            bar()
    if verbose: verbose_line_updater("Signals processed. Returning values.")
        
    # return the list.
    return repeat_id_remover(signal_id_list)



def query_speed_limits(geometry, limit_data_path, key="SPEED_LIM", epsg = 4326, margin = 1, last_known_limit= 20*0.00044704, verbose=False):
    '''
    query_speed_limits takes a geometric iterable of shapely points,
    a path to speed limit geodata, and returns the corresponding speed limits at
    each point.
    
    Params:
    geometry - a series of shapely points of lat and lon
    limit_data_path - path to shapefile of speed limit data as str
    key - column identifier of a speed limit or other value as str, assumed to be an int in units of mph
    margin - distance a point can be from the street to qualify in meters as int
    last_known_limit - the speed limit int value that the bus will assumedly start at or above. default 20 mph
    verbose - boolean to enable verbosity. Default False.
    
    Returns:
    list of speed limits corresponding to the geometry.
    
    Notes:
    11/13/2024 - This can be adjusted such that instead of the first value, it's the true closest.
                 Also can probably be combined with the other queries in such a way that it reduces individual
                 methods. FIXED -- 11/19/2024
    11/18/2024 - Currently Broken. - FIXED -- 11/19/2024
    '''
    # DEPENDING ON THE FILE, THIS CAN TAKE A WHILE.
    
    # get the bounding box of the route.
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))

    # get the speed limit data from the shapefile
    if verbose: verbose_line_updater("Loading shapefile & setting EPSG: {}".format(limit_data_path))
    limits = pio.read_dataframe(limit_data_path).to_crs(epsg=epsg)
    if verbose: verbose_line_updater("Shapefile loaded with EPSG:{}.".format(epsg))
    
    # get the streets that are within the bounding box
    if verbose: verbose_line_updater("Filtering streets...")
    bound_streets = limits[limits.apply(lambda x: shape_bounds.contains(x.geometry), axis=1).reset_index(drop=True) == True][['geometry', key]]
    if verbose: verbose_line_updater("Streets filtered.")

    # set up list.
    limit_list = []
    bound_streets['dist'] = 0
    # loop through the geometry:
    if verbose: verbose_line_updater("Processing speed limits...")
    with alive_bar(len(geometry)) as bar:
        for point in list(geometry):
            # get the speed limit at the point through checking that the distance is within the margin. Use first value that appears.
            bound_streets['dist'] = bound_streets.geometry.apply(lambda x: x.distance(point)*1000)

            # get the closest street
            closest = bound_streets[bound_streets['dist'] == bound_streets['dist'].min()].reset_index(drop=True)
            
            limit=0
            
            # if the distance is less than the margin, get the key value
            if closest['dist'][0] < margin:
                limit=closest[key][0]
            # append the limit to the list.
            limit_list.append(limit*0.00044704)
            bar()
    
    # while there is a speed limit of zero, swap them out for the next nearest speed limit that isn't zero
    
    while (0 in limit_list):
        remaining = limit_list.count(0)
        for i in range(len(limit_list)):
            if verbose: verbose_line_updater("Filling gaps: {}".format(remaining), reset=True)
            
            if limit_list[i] <= 0:
                limit_list[i] = last_known_limit
            else:
                last_known_limit = limit_list[i]
    if verbose: verbose_line_updater("Speed limits processed. Returning values.")


    # return the limit list, converted to km/s
    return limit_list


def query_bearings(geometry, bearing_type = "Angle", verbose=False):
    '''
    query_bearings takes a geometric iterable of shapely points,
    and returns the corresponding bearing when travelling form one point to the next.
    
    Params:
    geometry - a series of shapely points of lat and lon
    bearing_type - determine the bearing output.
                    Options: Angle - return the angle in degrees
                             Compass - return the angle as a directional bearing (eg: N, E, S, W)
    verbose - boolean to enable verbosity. Default False.
    
    Returns:
    list of bearings when travelling to subsequent point. 
    '''
    
    # get the current geometry series
    if verbose: verbose_line_updater("Loading geometry.")
    currents = geometry
    
    # shift the geometry series up by one, which will be the next ponts.
    if verbose: verbose_line_updater("Shifting geometry.")
    nexts = geometry.shift(-1)

    # combine the two into a dataframe.
    if verbose: verbose_line_updater("Combining geometry.")
    data = pd.concat([currents, nexts], axis=1, keys=['current', 'next']).reset_index(drop=True)

    # get the initial point and make sure it's the first value.
    final_point = list(data['current'])[-1]

    # set the nan value in the beginning of the 'nexts' to the initial, representing no distance traveled.
    data.iloc[-1, data.columns.get_loc("next")] = final_point

    # apply the haversine formula to the data, which will yeild a series of changes in distance.
    if verbose: verbose_line_updater("Calculating bearings.")
    bearings = data.apply(lambda x: point_bearing(x.current.y, x.current.x, x.next.y, x.next.x, bearing_type), axis=1)
    if verbose: verbose_line_updater("Bearings processed. Returning values.")
    
    return bearings

    
    

### ---- ELEVATION TOOLS ---- ###

def get_rasterfiles(dir_path):
    '''
    get_rasterfiles takes a path to a directory and returns a series of paths to any .tif files
    contained therin.
    
    Params:
    dir_path - path to directory, as str. 
    
    Returns:
    pandas series of path strings to each .tif file.
    '''
    
    # create a series using the list of items in the specified directory
    rasterfiles_raw = pd.Series(os.listdir(dir_path))
    
    # filter the list to only be those that contain tif. 
    rasterfiles_raw = rasterfiles_raw[rasterfiles_raw.apply(lambda x: ('tif' in (x.split("."))[-1]))].reset_index(drop='true')
    rasterfiles = rasterfiles_raw[rasterfiles_raw.apply(lambda x: ('.tif_' not in x))].reset_index(drop='true')
    # convert the paths to have the full path rather than just the filenames.
    rasterfiles = rasterfiles.apply(lambda x: "{}{}".format(dir_path, x))
    
    # return the series.
    return rasterfiles


def reproject_rasterfiles(filepath_sequence, target_crs='EPSG:4326', verbose=False):
    '''
    reproject_rasterfiles takes an iterable of geotiff files, and
    then reprojects them as a target crs and saves them in the same
    directory as the original.
    
    Params:
    filepath_sequence - an iterable of filepaths to geotiff files.
    target_crs - the projection system the rasterfiles are to be re-projected to.
    verbose - specify if the output is verbose. Boolean default False.
    
    Returns:
    iterable sequence of the newly reprojected files.
    
    Notes:
    11/14/2024 - This should be made to have the ability to override if
                 requested. For now, leaving it as-is is fine.
    '''
    
    # adjust the reprojection so it can be saved properly:
    crs_split = target_crs.split(":")
    crs_string = crs_split[0] + crs_split[1]
    
    if verbose: verbose_line_updater("Reprojecting raster files to: {}".format(target_crs))
    with alive_bar(len(filepath_sequence)) as bar:
        # Loop through each file in the iterable
        for path in list(filepath_sequence):

            #Check if the reprojection is already there.
            if (str(path.split(".tif")[1]) != '_{}'.format(crs_string) and (os.path.isfile(path + '_{}.tif'.format(crs_string)) == False)):

                # open the current rasterfile.
                with rasterio.open(path) as src:

                    # use rasterIO to calculate the transform of the file to the desired crs
                    transform, width, height = rasterio.warp.calculate_default_transform(src.crs,
                                                                                         target_crs,
                                                                                         src.width,
                                                                                         src.height,
                                                                                         *src.bounds)
                    # copy the kwargs from the geotiff
                    kwargs = src.meta.copy()

                    # update the kwargs with the target CRS
                    kwargs.update({
                        'crs': target_crs,
                        'transform': transform,
                        'width': width,
                        'height': height
                    })

                    # Create a reprojection file using the target crs
                    with rasterio.open(path+'_{}.tif'.format(crs_string), 'w', **kwargs) as dst:

                        # loop through the projections in src
                        for i in range(1, src.count + 1):

                            # reproject the values in the src.
                            rasterio.warp.reproject(source=rasterio.band(src, i),
                                                    destination=rasterio.band(dst, i),
                                                    src_transform=src.transform,
                                                    src_crs=src.crs,
                                                    dst_transform=transform,
                                                    dst_crs=target_crs,
                                                    resampling=rasterio.warp.Resampling.nearest)
            bar()
    if verbose: verbose_line_updater("Raster files succesfully reprojected.")
    
    return pd.Series(filepath_sequence).apply(lambda x: x + "_{}.tif".format(crs_string)).astype(str)


def query_elevation_series(geometry, gtiff_dir_ser, verbose=False):
    '''
    query_elevation_series takes a series of geometric shapely points, and an iterable of
    geotiff filepaths, and generates a pandas series of elevations from each point in
    the geometry.
    
    Params: 
    geometry - a series of shapely points of lat and lon
    gtiff_dir_ser - an iterable of geometric shapely points.
                    Must have an identical projection to the raster files.
    verbose - boolean parameter to specify verbosity. Defualt False.
    
    Returns:
    pandas series of elevations corresponding to the geometry iterable.
    '''
    
    # Generate the empty dictionary
    data_dict = {}
    
    # make sure the geometry is a series.
    geometry = pd.Series(geometry)
    
    if verbose: verbose_line_updater("Converting geometry...")
    # Convert the geometry to individual lat/lon tuples
    route_ll = (geometry.apply(lambda x: (x.x, x.y))).reset_index(drop=True)
    if verbose: verbose_line_updater("Geometry converted.")
    
    if verbose: verbose_line_updater("Querying elevation data...")
    # loop through each geotiff file in the series passed
    with alive_bar(len(gtiff_dir_ser)) as bar:
        
        for path in list(gtiff_dir_ser):

            # Use rasterio to open the file
            with rasterio.open(path) as img:

                # sample the elevations using rasterio and the list of lats and lons.
                route = pd.Series(list(rasterio.sample.sample_gen(img, list(route_ll), masked=True)))

                # remove any masked values
                filtered = route[route.apply(lambda x: str(type(x[0]))) != "<class 'numpy.ma.core.MaskedConstant'>"]
                
                # loop through the index of the filtered values
                for index in filtered.index:

                    # add the filtered value to the data dictionary at the index corresponding to
                    # its respective point.
                    data_dict[index] = filtered[index][0]
                bar()
                
    if verbose: verbose_line_updater("Elevation succesfully queried. Returning values.")
                
    # convert to series
    data_series = pd.Series(list(data_dict.values()), index=list(data_dict.keys()))
    
    # sort the data by index
    data_sorted = data_series.sort_index()
    
    # convert to km
    data_valued = data_sorted.apply(float).values/1000
    
    # return
    return data_valued




def smooth_elevation(elev_series, lg:int = 43, deg:int = 3):

        # Apply the savgol filter to the points
        y_new = scipy.signal.savgol_filter(elev_series,
                   lg,
                   deg,
                   axis = 0)

        return y_new


### ---- Processing Tools ---- ###


def combine_lidar_data(dsm, dx, max_grade=7.5):
    '''
    combine_lidar_data takes dsm elevation series, as well as the change in distance
    between the points, and generates a filtered elevation that provides cleaner road elevation.
    
    Params:
    dsm - digital surface model elevation iterable.
    dtm - digital terrain model elevation iterable.
    dx - iterable change in distance, in kilometers.
    max_grade - maximum possible grade the road can achieve. Default int=7.5.
    
    Returns:
    iterable of filtered elevations, in km.
    '''
    
    # combine the dsm and dx into a dataframe
    data = pd.concat([pd.Series(dsm), pd.Series(dx)], axis=1, keys=['dsm', 'dx']).reset_index(drop=True)
    
    # calculate the rolling dsm median
    data['rolling_dsm'] = (data['dsm'].rolling(7).median()).shift(-3)
    
    # calculate the difference in elevation
    data['dsm_dy'] = data['dsm'].diff()

    # get the filtered elevation. 
    filtered = data.apply(lambda x: (int(abs(x.dsm_dy/x.dx*100)>max_grade) * x.rolling_dsm) + (int(abs(x.dsm_dy/x.dx*100)<=max_grade)*x.dsm), axis=1)

    # return the values.
    return filtered.values


def calculate_grades(dx, elevations, clip = True , max_grade=7.5):
    '''
    calculate_grades takes the distance between points,
    as well as the elevations at each point, and returns the grade at each point.
    
    Params:
    dx - an iterable of distances between each point
    elevations - an iterable of elevations corresponding to each point.
    clip - boolean to determine to clip the grades or not. Default True.
    max_grade - an int representing the maximum grade a point can be without being clipped to. Default 7.5
    
    Returns:
    iterable of the slope grade at each point. 
    '''
    
    # calculate the grade percentage using the elevation difference and distance difference,
    grades = (pd.Series(elevations).diff()/pd.Series(dx)*100)
    
    # ensure only finite values
    grades = grades[np.isfinite(grades)]
    
    # check to clip the values
    if clip:
        grades = grades.clip(max_grade, -max_grade)
    else:
        grades=grades
    
    grades = list(grades)
    grades.append(grades[-1])
    # return the grade.
    return grades


def interpolate_geometry_series(geometry, max_distance=1):
    '''
    interpolate_geometry_series takes an iterable of geospatial shapely points,
    and then generates an interpolated series in the event that
    the maximum distance (in meters) is exceeded.
    
    Params:
    geometry - a series of shapely points of lat and lon
    max_distance - int representing maximum distance the points can be
                   without interpolation
                   
    Returns:
    iterable of lists of shapely points that have been interpolated.
    '''
    
    currents = pd.Series(geometry)
    
    # shift the geometry series up by one, which will be the next ponts.
    nexts = geometry.shift(-1)

    # combine the two into a dataframe.
    data = pd.concat([currents, nexts], axis=1, keys=['current', 'next']).reset_index(drop=True)
    
    # get the initial point and make sure it's the final value.
    final_point = list(data['current'])[-1]

    # set the nan value in the beginning of the 'nexts' to the initial, representing no distance traveled.
    data.iloc[-1, data.columns.get_loc("next")] = final_point
    # interpolate the points
    
    grace = data.apply(lambda x: interpolate_points(x.current, x.next, max_dist=max_distance), axis=1)
    
    # explode the data so that its a cohesive series
    sploded = grace.explode().reset_index(drop=True)
    
    # return an iterable of the interpolated points at each index
    return gpd.GeoSeries(sploded)


def interpolate_distance_traveled(dx, interpolated_iterable):
    
    data = pd.concat([pd.Series(dx), interpolated_iterable.apply(len)], axis=1, keys=['dx', 'i_i']).reset_index(drop=True)
    
    interp = data.apply(lambda x: [x.dx/x.i_i]*int(x.i_i), axis=1)
    
    return interp


        
        
    
### ---- Class, Saving, and Inter-operability tools ---- ###
class Route:
    '''
    Route class is is used to store route information. 
    
    Params:
    geometry - iterable of shapely geometry points
    elevation - iterable of elevation data
    limits - iterable of speed limit at each point. (optional, defaults to 25mph in km/s)
    stops - iterable of stop flags at each point. (optional, defaults to 10 evenly spaced stops.)
    signals - iterable of signal flags at each point. (optional, defaults to 8 evenly spaced signals.)
    signs - iterable of sign flags at each point. (optional, defaults to 3 evenly spaced stop signs.)
    
    Methods: 
    save_to_json() - saves the stored data to a json file. can be loaded again with Geography Tools' load_from_json.
    
    '''
    def __init__(self,
                 geometry,
                 elevation,
                 limits = None,
                 stops = None,
                 signals = None,
                 signs = None):
        
        # set up the params list
        self._params = [len(geometry),
                        not (limits is None),
                        not (stops is None),
                        not (signals is None),
                        not (signs is None)]
        
        # initialize geometry and elevation
        self.geometry = list(geometry)
        self.elevation = list(elevation)
        
        # get length of passed geometry
        geolen = len(geometry)
        
        # if limits is empty, use 25mph (in km/s) as default.
        if limits == None: self.limits = self._intersperse_list(25/2236.936, geolen)
        else: self.limits = list(limits)
        
        # if stops is empty, populate with -1, then intersperse with 10 stops (number of timepoints per route)
        if stops == None: 
            self.stops = self._intersperse_list(-1, geolen, True, 10)
            self.stops[-1] = 0
        else: self.stops = list(stops)
        
        # if signals is empty, intersperse with 10 signals.
        if signals == None: self.signals = self._intersperse_list(-1, geolen, 0, 8)
        else: self.signals = list(signals)
        
         # if signs is empty, intersperse with 3 signs.
        if signs == None: self.signs = self._intersperse_list(-1, geolen, 0, 3)
        else: self.signs = list(signs)
        
        # calculate the bearings, distance changes, and grades.
        self.bearings = query_bearings(pd.Series(self.geometry), bearing_type="Compass")
        self.dx = list(query_distance_traveled(pd.Series(self.geometry)))
        self.dz = query_elevation_changes(self.elevation)
        self.d_X = list(np.sqrt(np.asarray(self.dz)**2 + np.asarray(self.dx)**2))
        self.cum_d_X = list(pd.Series(self.d_X).cumsum())
        self.grades = calculate_grades(self.dx, elevation, max_grade=7.5)
        
        
        # return None
        return None
    
    
    def __str__(self):
        '''
        String method.
        '''
        return "Route(l:{},{},{},{},{})".format(self._params[0],
                                                self._params[1],
                                                self._params[2],
                                                self._params[3],
                                                self._params[4])
    
    
    def _intersperse_list(self, default_value, list_size, intersperse_val=-1, intersperse_num=0):
        '''
        intersperse_list() is used to generate a list of point-based flags of a given value,
        with the option to evenly intersperse those values with an alternate flag. 
        
        Params:
        default_value - the value that the list will be populated with initially.
        list_size - size of the list to be made, as an int.
        intersperse_val - value to be interspersed, default of -1.
        intersperse_num - number of ocurrences of the interspersed value, as int. Default 0.
        
        Returns:
        the new list.
        '''
        
        # Generate a new list of target size filled with target value.
        new_list = [default_value]*list_size
        
        # loop through and intersperse the desired value.
        for i in range(intersperse_num):
            new_list[int(i*(list_size/intersperse_num))] = intersperse_val
        
        # return the new list. 
        return new_list
    
    
    def save_to_json(self, path):
        '''
        save_to_json() is used to encode and save the Route data to a json file.
        
        Params:
        path - path and filename to be saved to, as str. Should end with '.json'.
        
        Returns:
        path to the saved json file.
        '''
        
        # generate data dictionary and compress the data where possible.
        data_dict = {'ge':encode_geometry(self.geometry),
                     'el':self.elevation,
                     'li':encode_series(self.limits),
                     'st':encode_series(self.stops),
                     'si':encode_series(self.signals),
                     'sn':encode_series(self.signs)}
        
        # Open and save to the json. encode in utf-8 for funsies.
        with open('{}'.format(path), 'w', encoding='utf-8') as f:
            
            # save the data.
            json.dump(data_dict, f, ensure_ascii=False, indent=4)
        
        # return the path.
        return path
    
    
    def save_as_routemap(self, path):
        '''
        save_as_routemap save the data as a classic routemap object's typical output.
        
        Params:
        path - path as str to save location. ends with .csv.
        
        Returns:
        path to CSV of saved routemap.
        
        Note:
        12/2/2024 - This is currently using a stand-in for ridership and ultimately should be PHASED OUT.
                    It is currently here so that it allows for backwards compatability with existing code.
        '''
        def ridership_generator(stops, mean_riders_per_trip=3.5, seed = 42):
            '''
            ridership generator takes a sequence of stop booleans,
            and the mean ridership for that trip, 
            and generates a sequence of the current riders at any given point.

            Params:
            stops - iterable of stop booleans
            mean_riders_per_trip - an int representing the mean riders per bus trip. default 3.5.
            seed - random seed for generation. Default is 42.

            Returns:
            list of ongoing ridership for a given point in the stop sequence.

            Notes:
            11/20/2024 - code breaks when mean riders requires more than 1 boarder per stop.
                         This should be remedied through the random generation system.
            '''
            random.seed(seed)
            # set up a list where the initial current riders is always 0
            current_riders = [0]

            print(stops.count(True))

            # calculate the maximum number of boarders
            range_max = int(mean_riders_per_trip/stops.count(True)*100)+1

            # loop through each stop
            for i in range(len(stops)):

                # if the stop is valid, generate the riders on and the riders off.
                if stops[i] == True:
                    riders_on = int(random.randrange(0, range_max+1) == range_max)
                    riders_off = random.randrange(0, current_riders[i]+1)

                    # update the rider list.
                    current_riders.append(current_riders[i] + riders_on - riders_off)
                else:
                    # update the riders list with the ongoing number of riders.
                    current_riders.append(current_riders[i])

            # remove the initial value.
            current_riders = current_riders[1:]

            # return the list.
            return current_riders

        
        stand_in_map = pd.concat([pd.Series(self.geometry),
                                  pd.Series(self.elevation),
                                  pd.Series(self.limits), # Convert to km/s
                                  pd.Series(self.stops).apply(lambda x: (x != -1)),
                                  pd.Series(self.signals).apply(lambda x: (x != -1)),
                                  pd.Series(self.dx)],
                                  axis=1, 
                                  keys=['geometry',
                                        'elevation[km]',
                                        'speed_limit[km/s]',
                                        'is_stop',
                                        'is_signal',
                                        'point_distances[km]']).reset_index(drop=True)
        stand_in_map['cumulative_distance[km]'] = stand_in_map['point_distances[km]'].cumsum()
        stand_in_map['latitude'] = stand_in_map['geometry'].apply(lambda x: x.x)
        stand_in_map['longitude'] = stand_in_map['geometry'].apply(lambda x: x.y)
        stand_in_map['ridership changes'] = pd.Series(ridership_generator(list(pd.Series(self.stops).apply(lambda x: bool(x+1))))).diff()
        stand_in_map.to_csv('{}'.format(path))
        
        return path
    
    
    def convert_to_RouteMap(self, path):
        '''
        convert_to_RouteMap() takes in a path for a saved RouteMap CSV (existing or not). If it
        exists, it loads it to a RouteMap object. If it doesn't, it converts the current Route
        object to a RouteMap object and saves at the given path. Used for Inter-operability with
        old code.
        
        Params:
        path - path to saved RouteMap csv file.
        
        Returns:
        RouteMap object.
        '''
        route_map = None
        if os.path.isfile(path):
            stand_in_map = pd.read_csv(path)
            route_map = rm.RouteMap(stand_in_map['geometry'], stand_in_map['elevation[km]'])
            route_map = route_map.load_from_dataframe(stand_in_map)
        else:
            self.save_as_routemap(path)
            stand_in_map = pd.read_csv(path)
            route_map = rm.RouteMap(stand_in_map['geometry'], stand_in_map['elevation[km]'])
            route_map = route_map.load_from_dataframe(stand_in_map)
        return route_map
    
    
    def query_point(self, index):
        '''
        query_point() takes an index of a point and returns the data for the corresponding point.
        
        Params:
        index - index of the route which data you would like to query.
        '''
        
        data = {'elevation':self.elevation[index],
                'limit':self.limits[index],
                'stop':self.stops[index],
                'signal':self.signals[index],
                'sign':self.signs[index],
                'bearing':self.bearings[index],
                'grade':self.grades[index],
                'dx_to_next':self.dx[index], # geodesic distance
                'dz_to_next':self.dz[index],
                'travel_dx':self.d_X[index],
                'sum_travel_dx':self.cum_d_X[index]}
        
        return data
    
    
    def to_gdf(self):
        '''
        to_gdf() takes the constituent information in geography_tools and converts it to a geodataframe.
        
        Params:
        N/A
        
        Returns:
        geodataframe of all information contained by a Route object.
        '''
        data = {'geometry':self.geometry,
                'elevation':self.elevation,
                'limit':self.limits,
                'stop':self.stops,
                'signal':self.signals,
                'sign':self.signs,
                'bearing':self.bearings,
                'grade':self.grades,
                'dx_to_next':self.dx, # geodesic distance
                'dz_to_next':self.dz,
                'travel_dx':self.d_X,
                'sum_travel_dx':self.cum_d_X}
        
        return gpd.GeoDataFrame(data)
        
            

def encode_series(iterable):
    '''
    encode_series() takes an iterable and compresses it down so that 
    repeating values are stored as the value and number of occurrences.
    Warning: Value typings will be lost.
    
    Params:
    iterable - an iterable list of values.
    
    Returns: 
    a string representation of the compressed information. 
    '''
    
    # get the first item of the iterable
    current_item = iterable[0]
    
    # set the current item count to negative to avoid indexing errors.
    current_item_count = -1
    
    # create an empty string.
    data_stream = ""
    
    # loop through each position in the iterable
    for i in range(len(iterable)+1):
        
        # if the index is at the end of the iterable,
        if i > len(iterable)-1:
            
            # add by 1 and subsequently add the current item and count to the string.
            current_item_count += 1
            data_stream += ('{},{}'.format(current_item_count, current_item))
        
        # otherwise,
        else:
            
            # the selected item is the current position in the iterable
            item = iterable[i]
            
            # and if the items are the same,
            if item == current_item:
                
                # increment the count by 1.
                current_item_count += 1
                
            # otherwise,
            else:
                
                # increment the count and add it to the stream.
                current_item_count += 1
                data_stream += ('{},{},'.format(current_item_count, current_item))
                
                # reset the count and set the new item.
                current_item_count = 0
                current_item = item

    # return the encoded representation.
    return data_stream


def decode_series(encoded_string):
    '''
    decode_series() takes in a string that was encoded using encode_series(), 
    and returns it to a decoded iterable. 
    
    Params:
    encoded_string - string that was encoded using encode_series()
    
    Returns:
    an iterable of the decompressed information. 
    '''
    
    # split the string by comma.
    spliterable = encoded_string.split(",")
    
    # create an empty list to contain the data.
    joint_iterable = []
    
    # for every two items, join them into a tuple and append it to the list.
    for item1, item2 in zip(spliterable[::2], spliterable[1::2]):
        joint_iterable.append((item1, item2))
        
    # create an empty list for the decoded information.
    decoded_iterable = []
    
    # loop through each item in the joint iterable
    for item in joint_iterable:
        
        # repeatedly add the encoded item for the number of occurences it has.
        for i in range(int(item[0])):
            decoded_iterable.append(item[1])
            
    # return the list. 
    return decoded_iterable


def encode_geometry(geo):
    '''
    encode_geometry() takes in an iterable of shapely points,
    and then converts it to a list of alternating x and y values.
    
    Params:
    geo - iterable of shapely points.
    
    Returns:
    list of alternating X and Y values in the same order as geo.
    '''
    
    # create an empty list
    clist = []
    
    # loop through the geometry and get the X and Y values,
    # and extend the list with them.
    for item in pd.Series(geo).apply(lambda x: [x.x, x.y]):
        clist.extend(item)
        
    # return the list.
    return clist


def decode_geometry(enc_geo):
    '''
    decode_geometry() takes in compressed geometry from encode_geometry(),
    and converts it to an iterable of shapely points.
    
    Params:
    enc_geo - encoded geometry list.
    
    Returns:
    list of shapely geometry.
    '''
    
    # empty list to store data
    point_list = []
    
    # loop through every 2
    for i in range(0, len(enc_geo), 2):
        # if even index, it's x, odd is y.
        # append the point.
        point_list.append(shapely.Point([enc_geo[i], enc_geo[i+1]]))
    
    # return the list.
    return point_list
            

def load_from_json(path):
    '''
    load_from_json() takes in a path to a json file from an exported Route object,
    and returns a path object as close as possible to what was saved.
    Waring: Some data types may be altered or adjusted. Speed
            limits are expected to be floats, stops, signals,
            signs, are expected to be ints.
    
    Params:
    path - path to json file.
    
    Returns: a Route object.
    '''
    
    # create empty data dictionary.
    data = {}
    
    # open the file.
    with open('{}'.format(path), 'r') as file:
        
        # load the data.
        data = json.load(file)
    
    # read the geometry data
    geo = decode_geometry(data['ge'])
    
    # read the elevation data
    el = data['el']
    
    # read the limit, stop, signal, and sign data
    lim = list(pd.Series(decode_series(data['li'])).apply(float))
    stp = list(pd.Series(decode_series(data['st'])).apply(int))
    signals = list(pd.Series(decode_series(data['si'])).apply(int))
    signs = list(pd.Series(decode_series(data['sn'])).apply(int))
    
    # generate a new Route object from the provided data. 
    loaded_route = Route(geo, el, lim, stp, signals, signs)

    return loaded_route
    
    
    