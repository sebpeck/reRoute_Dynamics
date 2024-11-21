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
    
        
    # return the stop id list with repeats removed.
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

    # get the initial point and make sure it's the last value.
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


def calculate_grades(dx, elevations, clip = False , max_grade=7.5):
    '''
    calculate_grades takes the distance between points,
    as well as the elevations at each point, and returns the grade at each point.
    
    Params:
    dx - an iterable of distances between each point
    elevations - an iterable of elevations corresponding to each point.
    clip - boolean to determine to clip the grades or not. Default False.
    max_grade - an int representing the maximum grade a point can be without being clipped to.
    
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


        
        
    
    