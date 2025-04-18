'''
KC_Query_Tools.py
KC_Query_Tools is used for providing route information from King County Metro's dataset.
'''
import sys
sys.path.append('../src/')
import pandas as pd
import geopandas as gpd
import shapely
import numpy as np
import pyogrio as pio
from reRoute_Dynamics import Geography_Tools as gt
import os
import multiprocessing
import datetime as dt
import inspect


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
                                                                "KC_Query_Tools",
                                                                str(inspect.stack()[1][3]),
                                                                message)
    # check if the carrage should be reset.
    if reset:
        print(data_string, end='\r')
    else:
        print(data_string)
        
    # return the carrage.
    return data_string


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
    

def _query_shape_geometry(shape_num, shape_table_dir, fwd_flag = True):
    '''
    query_shape_geometry takes in a shape number, shape lookup directory, and returns the 
    corresponding geometry in lat, lon
    
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
        return gdf['geometry'].apply(lambda x: shapely.ops.transform(flip, x))
    else:
        return gdf.iloc[::-1]['geometry'].apply(lambda x: shapely.ops.transform(flip, x))
    
    
    
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


def _flip(x, y):
    """Flips the x and y coordinate values"""
    return y, x


def render_kc_route_file(fname, saved_route_list, route_data_dir, dtm_raster_path, render_params):
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
    render_params - the parameters as a touple to help define interpolation and smoothness degree
    
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
            geo = geo.apply(lambda x: shapely.ops.transform(_flip, x))
            
            dx = gt.query_distance_traveled(geo, verbose=verb)

            # use the geometry to interpolate points
            i_geo = gt.interpolate_geometry_series(geo, render_params[0])

            # Use the interpolated geometry to get interpolated distances
            i_dx = gt.query_distance_traveled(i_geo, verbose=verb)

            # Get the smoothness based on the interpolation ratio
            smoothness=int(np.ceil((len(i_dx)/len(dx))**2))+1
            
            if len(render_params) > 2:
                smoothness=render_params[2]

            # get the DTM elevation and sm0oth it.
            idtm_elevation = gt.query_elevation_series(i_geo, reprojected_DTMs, verbose=verb)*0.3048
            smellev = gt.smooth_elevation(idtm_elevation,smoothness, render_params[1])

            # get the speed limits for the data
            limits = gt.query_speed_limits(i_geo, route_data_dir + "Seattle_Streets/bound_KCM_Streets.shp", verbose=verb)

            # Get the stops for the data
            stops = gt.query_stops(i_geo, route_data_dir + "stops.txt", verbose=verb)

            # validate the stops
            stops = check_valid_stops_by_shape(stops, shape_id, route_data_dir + "trips.txt", route_data_dir + "stop_times.txt")

            # make the ends be zero - Why was this the case??
            stops[0] = 0
            stops[-1] = 0

            # Get the signal data
            signals = gt.query_signals(i_geo, route_data_dir +"Signals/bound_KCM_signals.shp", verbose=verb)

            # Generate a new test route
            route = gt.Route(i_geo, smellev, limits=limits, stops=stops, signals=signals, signs = [-1]*len(i_geo), smooth_grades=True )
            route.save_to_json(fname)
            print("Succesfully Rendered {}".format(fname))
        return fname + ""
    except Exception as e:
        import traceback
        print("Failed Rendering {}, Error: {}".format(fname, e))
        traceback.print_exc()
        return fname + " -- 1 "
    
    
def batch_render_kc_routes(route_options, route_data_dir, elevation_raster_path, route_savepath, skip_unrendered=False, render_params=(10, 3, 25),batch_size = 5, verbose=False):
    '''
    batch_rended_kc_routes takes in a list of route shortnames, path to king county route data,
    a path to an elevation raster, and a path to write the saved .json files to, and then
    renders all the given route options, both inbound and outbound, for each possible shape of
    that shortname.
    
    Params:
    route_options: list of route shortnames, all as ints
    route_data_dir: path to directory, as str, that holds route data from king county (a la, routes.txt, trips,txt.)
    elevation_raster_path: path to directory containing all elevation rasters for the area being examined. 
    route_savepath: path to directory the json files will be saved to.
    skip_unrendered: boolean (default False) to skip rendering and just return the pre-existing valid json files.
    render_params: tuple containing the parameters for (interpolation distance, savgol degree) as ints. Default (10, 3)
    batch_size: (default 5) - int representing how many batches/cpu cores to use while rendering.
    verbose: boolean to enable verbosity.
    
    Returns:
    a list to hold the paths of all files that have been processed. 
    If the render failed, the path will have ' -- 1' added to the end.
    
    '''

    # generate the route list of saved routes that exist
    saved_route_list = os.listdir(route_savepath)
    for i in range(len(saved_route_list)):
        saved_route_list[i] = (route_savepath+saved_route_list[i])


    # empty list to hold the json filenames that need to be generated
    fname_list = []

    # empty list to hold valid files that have been generated
    data_list = [[]]


    # loop through each route option, and iterate through the shapes and directions
    # to hold the json filenames
    for selected_route in route_options:
        # try to load the selected route by shortname
        try:
            # Query the possible shapes of that route.
            route_id = query_route_id(selected_route, route_data_dir + "routes.txt")
            possible_shapes = query_possible_shapes(route_id, route_data_dir + "trips.txt")


            # loop through each possible shape
            for shape_id in possible_shapes:

                # loop through inbound and outbound
                for in_out in [0, 1]:

                    # generate the json filename based on shortname, shape, and direction
                    filename = route_savepath + 'rt{}_sh{}_d{}.json'.format(selected_route, shape_id, in_out)
                    # Check if the filename doesn't exist in the route list
                    if filename not in saved_route_list:

                        # if not, add the name to the list of json files that need to be generated
                        fname_list.append(filename)

                    else:
                        # if it does exist, add it to the first list of valid files in the data list.
                        data_list[0].append(filename)
                        if verbose: verbose_line_updater('Filename {} exists already. Skipping...'.format(filename))
        except Exception as e:
            if verbose: verbose_line_updater("Skipping route {}... Possible error: {}".format(selected_route, e))

    
    if not skip_unrendered:

        # batch the filenames based on the batch size
        batched_filenames = [fname_list[i:i + batch_size] for i in range(0, len(fname_list), batch_size)]


        # loop through each batch of filepaths
        for fname_batch in batched_filenames:

            # Get the batch file and associate it with its parameters by zipping to touples
            batch = list(zip(fname_batch,
                     [saved_route_list]*len(fname_batch),
                     [route_data_dir]*len(fname_batch),
                     [elevation_raster_path]*len(fname_batch),
                     [render_params]*len(fname_batch)))

            # use multiprocessing with a batch size for the pool
            with multiprocessing.Pool(batch_size) as pool:
                # use starmap to process each file's route map.
                data_list.append(pool.starmap(render_kc_route_file, batch))
    
    return data_list


def calculate_expected_ridership(ridership_data_path):
    '''
    calculate_expected_ridership takes a path to king county ridership data by route, 
    and then calculates the typical total number of unique riders for a given route,
    period, and direction. 
    
    Params:
    ridership_data_path: Path, as str, to ridership data csv
    
    Returns:
    dataframe containing the total unique riders for that given combination of route, period, and direction.
    '''
    
    # Load the CSV
    ridership_frame = pd.read_csv(ridership_data_path)

    # Group the data by route, period, and direction
    grouped = ridership_frame.groupby(['Route', 'Period', 'InOut'])
    formatted_rider_data = []
    
    # iterate through the grouped data
    for combo_group_name, combo_group in grouped:
        
        # group each combined group by the stop sequence
        sid_groups = combo_group.groupby('STOP_SEQ')
        
        # create a list to store the data
        data_list = []
        
        # Loop through eaach sequence:
        for stop_id in sid_groups.groups:
            
            # Get the group of sequences
            stop_group = sid_groups.get_group(stop_id)
            
            # get the mean of the average on, off, and load
            ave_on = stop_group['AveOn'].mean()
            #ave_off = stop_group['AveOff'].mean()
            #ave_Ld = stop_group['AveLd'].mean()
            
            # convert the information to a data dict
            data = {'STOP_SEQ':stop_id, 'diff':ave_on}
            
            # append the dict to the data list
            data_list.append(data)
            
        # generate a dataframe from the data list
        g_in = pd.DataFrame(data_list)
        
        # get the cumulative value of riders
        expected_unique_riders = g_in['diff'].cumsum().max()
        
        # Get the data to be exportde
        group_data = {'rt':combo_group_name[0],
                      'per':combo_group_name[1],
                      'io':combo_group_name[2],
                      'riders':expected_unique_riders}
        
        # append it to the formatted list
        formatted_rider_data.append(group_data)
        
    # generate a datagrame of the ridership data
    formatted_rider_data = pd.DataFrame(formatted_rider_data)
    
    #return said dataframe
    return formatted_rider_data

def get_ridership(shortname, period, rider_data):
    '''
    get_ridership takes a shortname, period, and ridership dataframe, 
    and provides a list of dictionaries, with each key being period,
    and the value being unique ridership for that period. 
    '''
    return list(rider_data.groupby(['rt', 'io']).get_group((int(shortname), period))[['per', 'riders']].apply(lambda x: {x['per']: x['riders']}, axis=1))


def bus_type_finder(val):
    '''
    bus_type_finder takes a king county bus ID number, 
    and returns the class of bus it belongs to as a string. 
    '''
    if 6813<=val<=6865:
        return 'DE60LF'
    elif 6000<=val<=6019:
        return 'DE60LFA'
    elif 7001<=val<=7199:
        return 'OBI6'
    elif 6866<=val<=6999 or val==6800 or 6020<=val<=6035 or 6040<=val<=6073 or 6075<=val<=6117:
        return 'DE60LFR'
    elif 3700<=val<=3759:
        return 'XDE35'
    elif 7200<=val<=7259:
        return 'XDE40'
    elif 4300<=val<=4409:
        return 'XT40'
    elif 4500<=val<=4563:
        return 'XT60'
    elif 6200<=val<=6269 or 8000<=val<=8084 or 8100<=val<=8299 or 6400<=val<=6412:
        return 'XDE60'
    elif 7300<=val<=7484:
        return 'BRT'
    elif 4700<=val<=4719:
        return 'XE40'
    elif 4800<=val<=4819:
        return 'XE60'
    else:
        return "NA"
