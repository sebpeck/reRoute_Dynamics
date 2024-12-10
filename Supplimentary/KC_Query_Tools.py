'''
KC_Query_Tools.py
KC_Query_Tools is used for providing route information from King County Metro's dataset.
'''

import pandas as pd
import geopandas as gpd
import shapely
import numpy as np
import pyogrio as pio


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

