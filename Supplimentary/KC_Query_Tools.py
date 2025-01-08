'''
KC_Query_Tools.py
KC_Query_Tools is used for providing route information from King County Metro's dataset.
'''

import pandas as pd
import geopandas as gpd
import shapely
import numpy as np
import pyogrio as pio
import Geography_Tools as gt


def query_route_id(shortname, route_data_path, shortname_key = 'route_short_name', idnum_key = 'route_id'):
    '''
    query_route_id takes in a route's shortname and path to the compiled route datafile,
    and returns the shortname's route_id as an int.
    
    Params:
    shortname - route shortname as str
    route_data_path - path to route data csv as str
    shortname_key - key for identifying shortname column in the data as str, default 'route_short_name'
    idnum_key - key for identifying route id column in the data as str, default 'route_id'
    
    Returns:
    route_id corresponding to shortname.
    '''
    
    # load the route id table
    route_id_table = pd.read_csv(route_data_path)
    
    # load the route list from the table
    route_list = list(route_id_table[shortname_key])
    
    # get the route id from the data by filtering to only the routes with the shortname and getting the first.
    route_id = (route_id_table[route_id_table[shortname_key] == str(shortname)].reset_index(drop=True))[idnum_key][0]
    
    # return the route id
    return route_id


def query_possible_shapes(route_id, shape_data_path, idnum_key = 'route_id', shapeid_key = 'shape_id'):
    '''
    query_routes takes a route id and the path to the shape data, 
    and returns a list of the possible shape ids for that route.
    
    Params:
    route_id - id corresponding to a route in the shape data
    shape_data_path - path to shape id data csv as str
    idnum_key - key for identifying route id column in data as str, default 'route_id'
    shapeid_key - key for identifying shape id column in data as str, defualt 'shape_id'
    
    Returns:
    numpy array of possible shape ids for given route id
    '''
    
    # load the data
    trip_table = pd.read_csv(shape_data_path)
    
    # filter the trip table to the trips the route takes
    route_trips = trip_table[trip_table[idnum_key] == route_id].reset_index(drop=True)
    
    # using the possible treps, get all the unique shape identifiers
    possible_shapes = route_trips[shapeid_key].unique()
    
    # return the array of possible shapes
    return possible_shapes


def query_shape(shape_num, shape_table_dir, fwd_flag = True):
    '''
    query_shape takes in a shape number, shape lookup directory, and returns the 
    corresponding geodataframe of that sequence. Distance is in km.
    
    Params:
    shape_num - number corresponding to a shape
    shape_table_dir - path to file containing shape geodata csv as str
    fwd_flag - boolean to report the shape in forward or reverse. Default is True (forward)
    
    Returns:
    geodataframe corresponding to the provided data. Distances in km. 
    '''
    
    # load the data
    shape_table = pd.read_csv(shape_table_dir)
    
    # get the shape corresponding to the shape id through filtering
    trip_shape = shape_table[shape_table['shape_id'] == shape_num]
    
    # convert the lat and lon to shapely points.
    trip_geom = trip_shape.apply(lambda x: shapely.Point(x.shape_pt_lon, x.shape_pt_lat), axis=1)
    
    # convert the baseline distance traveled to kilometers
    trip_dist = trip_shape['shape_dist_traveled']*0.0003048
    
    # convert the data to a geodataframe
    gdf = gpd.GeoDataFrame(data=trip_dist, geometry=trip_geom)
    
    # add the sequence indexer
    gdf['sequence'] = trip_shape['shape_pt_sequence']
    gdf = gdf.reset_index(drop=True)
    
    # Check if the series should be reversed or not before returning
    if fwd_flag:
        return gdf
    else:
        return gdf.iloc[::-1]
    
    
def check_valid_stops_by_shape(stop_sequence, shape_id, trip_table_dir, stop_times_dir):
    '''
    check_valid_stops_by_shape() takes an iterable of stop IDs,
    a shape id, and the directories for trip table and stop time data,
    and validates all stops in that sequence that belong to the passed shape.
    
    Params:
    stop_sequence - an iterable of stop ID ints.
    shape_id - an int representing a route shape.
    trip_table_dir - str representing path to the trip table.
    stop_times_dir - str representing path to stop times data.
    
    Returns:
    list of booleans corresponding to valid or invalid stops.
    '''
    
    # empty list to generate
    stop_bools = []
    
    # read the data
    trips = pd.read_csv(trip_table_dir)
    stop_times = pd.read_csv(stop_times_dir)
    
    # loop through each id in the passed list
    for stop_id in stop_sequence:   
        
        # filter the stop times to the ones with the current id
        s_stop_times = stop_times[stop_times['stop_id'] == stop_id]
        
        # empty list for possible shapes that use that stop
        shape_id_list = []
        
        # loop through the unique trip ids
        for trip_id in s_stop_times['trip_id'].unique():
            
            # add all unique shapes corresponding to that trip to the shape list
            shape_id_list.extend(trips[trips['trip_id'] == trip_id]['shape_id'].unique())
            
        # gives true if the shape ID is within the uniqe values pulled from the trips with that stop.
        is_valid_stop = shape_id in pd.Series(shape_id_list).unique()
        
        # append the stop boolean
        stop_bools.append(int(is_valid_stop)-1)
        
    # return the list.
    return stop_bools



def render_kc_route_file(fname, saved_route_list, route_data_dir, dtm_raster_path, trip_obj):
    '''
    render_kc_route_file is the pipeline for rendering king county specific routes.
    
    Params:
    fname - name of route json file to be exported. Should be in format of 'rt{}_sh{}_d{}.json'
            in order to be parsed for shape ID.
    saved_route_list - list of json files in directry to check wether the given fname has been saved yet.
    route_data_dir - path to directory containing the route data, which should contain the following:
                     'shapes.txt' -- a file containing the KCM shapes
                     'Seattle_Streets/KCM_Streets.shp' shapefile of all king county metro streets
                     'stops.txt' -- a file containing the KCM bus stops
                     'trips.txt' -- a file containing the KCM trip data
                     'stop_times.txt' -- a file containing the stop time KCM data
                     'Signals/KCM_signals.shp' -- a shapefile containing the signals for KCM area
    dtm_raster_path: path to raster elevation data, should be in same projection as all other data (EPSG:4326)
    trip_obj - the trip parameters object to help define interpolation and smoothness degree
    
    Returns:
    filename of object with a ' -- 0' or ' -- 1' appended, depending on if the render succeeded or not.
                    
    '''
    verb = False
    try:
        if fname not in saved_route_list:
            print("Rendering File {}".format(fname))
            shape_id = fname.split('_sh')[1]
            in_out = int(shape_id.split('_d')[1].split('.json')[0])
            shape_id = int(shape_id.split('_')[0])
            
            # Using geography_tools, process and reproject any DTMs. 
            DTM_rasters_list= gt.get_rasterfiles(dtm_raster_path)#+ "_reproject.tif"
            DTM_rasters_list = DTM_rasters_list[DTM_rasters_list.apply(lambda x: 'reproject' in x) == False].reset_index(drop=True)
            reprojected_DTMs = gt.reproject_rasterfiles(DTM_rasters_list, verbose=verb)
            # load the shape
            shape = query_shape(shape_id, route_data_dir + "shapes.txt", bool(in_out))

            # get the geometry and distances
            geo = shape['geometry']
            dx = gt.query_distance_traveled(geo, verbose=verb)

            # use the geometry to interpolate points
            i_geo = gt.interpolate_geometry_series(geo, trip_obj.d_interp)

            # Use the interpolated geometry to get interpolated distances
            i_dx = gt.query_distance_traveled(i_geo, verbose=verb)

            # Get the smoothness based on the interpolation ratio
            smoothness=int(np.ceil((len(i_dx)/len(dx))**2))+1

            # get the DTM elevation and smoth it.
            idtm_elevation = gt.query_elevation_series(i_geo, reprojected_DTMs, verbose=verb)
            smellev = gt.smooth_elevation(idtm_elevation,smoothness, trip_obj.deg)

            # get the speed limits for the data
            limits = gt.query_speed_limits(i_geo, route_data_dir + "Seattle_Streets/KCM_Streets.shp", verbose=verb)

            # Get the stops for the data
            stops = gt.query_stops(i_geo, route_data_dir + "stops.txt", verbose=verb)

            # validate the stops
            stops = check_valid_stops_by_shape(stops, shape_id, route_data_dir + "trips.txt", route_data_dir + "stop_times.txt")

            # make the ends be zero - Why was this the case??
            stops[0] = 0
            stops[-1] = 0

            # Get the signal data
            signals = gt.query_signals(i_geo, route_data_dir +"Signals/KCM_signals.shp", verbose=verb)

            # Generate a new test route
            route = gt.Route(i_geo, smellev, limits=limits, stops=stops, signals=signals, signs = [-1]*len(i_geo))
            route.save_to_json(fname)
            print("Succesfully Rendered {}".format(fname))
        return fname + ""
    except:
        print("Failed Rendering {}".format(fname))
        return fname + " -- 1"


