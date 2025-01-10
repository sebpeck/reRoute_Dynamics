'''
KC_Route_Testing.py
'''
import sys
sys.path.append('../src/reRoute_Dynamics_Core/')
import pandas as pd
import numpy as np
import datetime
import multiprocessing
import Trip_Simulator as ts
import KC_Query_Tools as kqt
import Object_Params as op
import Geography_Tools as gt
import ESS
import matplotlib.pyplot as plt


def prepare_trips(data_list,
                 formatted_rider_data,
                 trip):
    '''
    prepare_trips takes a list of valid route save paths, formatted ridership data,
    and then combines the two into a format that's easily iterable with an individual trip object for each 
    path, not including seeds.
    
    Params:
    data_list: list of str paths to json files as exported by render_kc_route_file or batch_render_files (from KC_Query_Tools.py) 
    formatted_rider_data: pandas dataframe containing data exported by calculate_expected_ridership (from KC_Query_Tools.py)
    trip: trip object from Object_Params
    
    Returns:
    dataframe of prepped trips.
    '''
    
    # Get all the saved routes that are valid by:
    # Making an empty list
    saved_routes = []

    # looping through the data list to verify validity of export
    for listed_routes in data_list:

        # Check if the entry has an error
        valid_routes = [ x for x in listed_routes if " -- 1" not in x ]
        valid_routes = list(pd.Series(valid_routes))

        # extend the saved routes with only the valid ones.
        saved_routes.extend(valid_routes)

    # generate a dataframe of the saved routes
    routes = pd.DataFrame(pd.Series(saved_routes))

    # Make the column 'path'
    routes.columns = ['path']

    # split the paths in such a way that you can extract the shortname, shape, and directory
    route_id_lists = routes['path'].apply(lambda x: x.split('/')[-1]).apply(lambda x: x.split('.json')[0]).apply(lambda x:x.split('_'))
    routes['route'] = route_id_lists.apply(lambda x: x[0][2:])
    routes['shape'] = route_id_lists.apply(lambda x: x[1][2:])
    routes['dir'] = route_id_lists.apply(lambda x: x[2][1:])

    # creat a list to hold the paths that are actually testable based on the info we have ridership for
    paths_list = []
    # Loop through the ridership information
    for index, row in formatted_rider_data.iterrows():

        # Get the ridership shortname, period, and direction
        route_name = row['rt']
        period = row['per']
        io = row['io']
        riders = row['riders']

        # Convert direction to a 0 or 1
        if io == 'I':
            io = 0
        elif io == 'O':
            io = 1

        # get the json files that match the route and direction.
        paths = routes[routes['route'] == str(route_name)]
        paths = paths[paths['dir'] == str(io)]
        paths = list(paths['path'])

        # for each file that matches the route and direction,
        for path in paths:
            # create a new entry of shortname, period, direction, ridership, and path.
            modified_trip = trip.copy()
            modified_trip.m_riders = riders
            
            # export the entry data
            entry = {'rt':route_name, 'per':period, 'io':io, 'trip':modified_trip, 'path':path}
            paths_list.append(entry)

    # convert the paths to be rendered to a dataframe
    testing_routes = pd.DataFrame(paths_list)
    return testing_routes


def run_trip(path, trip, bus, bus_ESS, export_figures = False, seed_list=np.arange(0,5,1)):
    '''
    run_trip takes a path to a route savefile json, 
    a trip object, bus object, and ESS object, 
    and runs the given trip, returning the miles/kwh results.
    
    Params:
    path: path, as str, to a json file of a route savefile as exported by KC_Query_Tools
    trip: trip object, as from Object_Params
    bus: bus object, as from Object_Params
    bus_ESS: ESS object, as from Object Params
    seed_list: list of seeds to be used. Trip will be run once per seed. Default [0, 1, 2, 3, 4]
    export_figures: default false. Determines if figures should be rendered for each. THIS CAN BE VERY MEMORY INTENSIVE.
    
    Return:
    list of lists of mi/kwh values, with each index of the exterior list corresponding with a single trip combo,
    and each index in the sublist correspinding to each seed. 
    IF EXPORT_FIGURES is TRUE, this will export as a tuple, with 0 being the aformentioned, and 1 being
    the same format, but with the plots instead. THIS DOES NOT WORK WITH MULTIPROCESSING.
    '''
    
    # empty list for results.
    results = []
    plt_list = []
    
    
    # Loop through each seed.
    for seed in seed_list:
        
        # replace the trip seed with the current seed. 
        trip.seed = int(seed)
        
        # set up a last_gdf variable
        last_gdf = None
        
        # load the route.
        test_route = gt.load_from_json(path)
        
        # set the last gdf.
        last_gdf = test_route.to_gdf()

        # simulate the trip to get the running data.
        running_data = ts.simulate_trip(test_route, trip = trip, bus = bus, ESS = bus_ESS)
        
        # Convert the running data to a dataframe.
        route_results = pd.DataFrame(running_data)

        # Use ESS to calculate the battery power from the required power. 
        route_results['BP'] = route_results['P'].apply(lambda x: ESS.calc_instance_power(x,
                                                                                         motor_eff=bus_ESS.Em,
                                                                                         invert_eff=bus_ESS.Ei,
                                                                                         aux_eff=bus_ESS.Ea,
                                                                                         aux_load=bus_ESS.P_aux,
                                                                                         regen_eff=bus_ESS.Er,
                                                                                         max_regen=bus_ESS.P_regen))
        
        
        
        
        #route_results[(route_results['BP']).apply(abs) < 10000]

        # Convert the distances to miles from meters
        mi = route_results['dx'].sum()/1609.344
        
        # convert the net energy to kwh.
        net_energy = (route_results['BP']/1000*route_results['dt']).sum()/3600
        
        # append mi/kwh to the results
        results.append(mi/net_energy)
        
        # if we're exporting figures,
        if export_figures:
            
            # generate the stop markers
            stop_markers = route_results.apply(lambda x: (sum(x['stop_clf'])>0), axis=1)
            
            # generate the figure
            fig, ax = plt.subplots(3,2, figsize = (17, 12), dpi=300, gridspec_kw={'width_ratios': [3, 1]})

            # plot the geodesic distance vs the battery power
            ax[0,0].plot(route_results['gdx'].cumsum(), route_results['BP'], color = '#b319b3')

            ax[0,0].set_zorder(1)
            ax[0,0].set_facecolor('none')
            ax[0,0].set_xlabel('Geodesic Distance [m]')
            ax[0,0].set_ylabel('ESS Power [Watts]')
            
            # plot geodesic distance vs elevation
            ax2_0_0 = ax[0,0].twinx()
            ax2_0_0.fill_between(route_results['gdx'].cumsum(), route_results['elevation'], color = '#BDBDBD')
            ax2_0_0.scatter(route_results['gdx'].cumsum(), stop_markers*route_results['elevation'], color='tab:red')
            ax2_0_0.set_ylim(0.01)
            ax2_0_0.set_zorder(0)
            ax2_0_0.set_ylabel('Elevation [km]')

            # plot the histogram of distance traveled at a given power. 
            bins = [-140,-120,-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
            ax[0, 1].hist(pd.Series(route_results['BP'])/1000, color='#b319b3', bins=bins, weights=route_results['dx'])
            ax[0,1].set_xlabel('ESS Power [kW]')
            ax[0,1].set_ylabel('Cumulative distance [m]')

        
            # plot the battery power vs time
            ax[1,0].plot(route_results['dt'].cumsum(), route_results['BP'], color = '#b319b3')

            ax[1,0].set_zorder(1)
            ax[1,0].set_facecolor('none')
            ax[1,0].set_xlabel('Elapsed Time [s]')
            ax[1,0].set_ylabel('ESS Power [Watts]')
            
            # Plot the time vs elevation
            ax2_1_0 = ax[1,0].twinx()
            ax2_1_0.fill_between(route_results['dt'].cumsum(), route_results['elevation'], color = '#BDBDBD')
            ax2_1_0.scatter(route_results['dt'].cumsum(), stop_markers*route_results['elevation'], color='tab:red')
            ax2_1_0.set_ylim(0.01)
            ax2_1_0.set_zorder(0)

            # plot the histogram of time spent at a given power
            bins = [-140,-120,-100, -80, -60, -40, -20, 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
            ax[1, 1].hist(pd.Series(route_results['BP'])/1000, color='#b319b3', bins=bins, weights=route_results['dt'])
            ax[1,1].set_xlabel('ESS Power [kW]')
            ax[1,1].set_ylabel('Cumulative Time [s]')

            # plot the histogram of time spent at a given grade.
            bins2 = [-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8]
            ax[2, 1].hist(pd.Series(route_results['grade']), color='#BDBDBD', bins=bins2, weights=route_results['dt'])
            ax[2,1].set_xlabel('grade (%)')
            ax[2,1].set_ylabel('Cumulative Time [s]')

            # plot the geodesic distance vs velocity. 
            ax[2, 0].plot(route_results['gdx'].cumsum(), route_results['v_f'])
            ax[2,0].scatter(route_results['gdx'].cumsum(), (stop_markers), color='tab:red')

            # append the plot to the list
            plt_list.append(plt)
            
    # if exporting figures, export the tuple
    if export_figures:
        return (route_results, plt_list)

    # otherwise, export just the results.
    else:
        return results


def hyperbaric_time_chamber(testing_routes,
                            bus,
                            bus_ess,
                            seed_list = np.arange(0,5,1),
                            batch_size = 5,
                            end_dex=0):
    '''
    hyperbaric_time_chamber takes a dataframe of valid testing routes, an ESS, 
    
    Params:
    testing_routes: a dataframe of 
    trip: trip object from Object_Params
    bus: bus object from Object_Params
    bus_ess: ESS object from Object_Params
    seed_list: seeds to iterate through, as a list. Default [0,1,2,3,4].
    batch_size: (default 5) number of items per batch, or CPU cores to be used when processing. 
    end_dex: stopping point, default uses all/value 0. Used for running small tests and debugging. 
    
    Returns:
    list containing the list of results for each combination of period and route shape.
    '''

    # create batches of route paths and trips
    batches_of_routes = [list(testing_routes['path'])[i:i + batch_size] for i in range(0, len(list(testing_routes['path'])), batch_size)][:end_dex]
    batches_of_trips = [list(testing_routes['trip'])[i:i + batch_size] for i in range(0, len(list(testing_routes['trip'])), batch_size)][:end_dex]
    
    if end_dex != 0:
        batches_of_routes = batches_of_routes[:end_dex]
        batches_of_trips = batches_of_trips[:end_dex]

    # create a lsit to hold results
    route_milage_data = []
    
    # for each route in the batches
    for i in range(len(batches_of_routes)):
        print("{}/{}, {}".format(i,len(batches_of_routes), datetime.datetime.now().time()))
        route_batch = batches_of_routes[i]
        trip_batch = batches_of_trips[i]
        full_batch = zip(route_batch, trip_batch, [bus]*len(route_batch), [bus_ess]*len(route_batch), [False]*len(route_batch), [seed_list]*len(route_batch))
        with multiprocessing.Pool(batch_size) as pool:
            route_milage_data.extend(pool.starmap(run_trip, full_batch))

    return route_milage_data

def run_tests(ridership_data_path,
              route_data_dir,
              elevation_raster_path,
              route_savepath,
              bus_path,
              ess_path,
              trip_path,
              batch_size=5,
              seed_list = np.arange(0, 5, 1),
              end_dex = 0):
    '''
    run_tests runs load tests, specifically formatted for King County's data. 
    
    Params:
    ridership_data_path: path to KC ridership data
    route_data_dir: path to KC route data directory
    elevation_raster_path: path to elevation raster data
    route_savepath: path to json savefiles for route geodata
    bus_path: path to saved bus object
    ess_path: path to saved ess object
    trip_path: path to saved trip object
    batch_size: default 5, int representing num of cores to use while multiprocessing
    seed_list: list of seeds to use when iterating, default [0,1,2,3,4]
    end_dex: default 0 - used for debugging by only running a limited number of the batches.
    
    Returns:
    dataframe containing results!
    
    '''
    # Load in the objects
    bus = op.load_bus_params(bus_path)
    ess = op.load_ESS_params(ess_path)
    trip = op.load_trip_params(trip_path)
    
    # extract ridership data for full factorial combinations
    formatted_rider_data = kqt.calculate_expected_ridership(ridership_data_path)

    # And also get all the routes with ridership data while converting the values to string for inter-operability.
    route_options = list(pd.Series(formatted_rider_data['rt'].unique()).apply(str))

    # render out the files. 
    data_list = kqt.batch_render_kc_routes(route_options,
                                           route_data_dir,
                                           elevation_raster_path,
                                           route_savepath,
                                           skip_unrendered=True,
                                           render_params=(trip.d_interp, trip.deg),
                                           batch_size=batch_size,
                                           verbose=True)
    
    prepped_trip_frame = prepare_trips(data_list, formatted_rider_data, trip)
    trip_results = hyperbaric_time_chamber(prepped_trip_frame,
                                            bus,
                                            ess,
                                            seed_list = seed_list,
                                            batch_size = batch_size,
                                            end_dex=end_dex)
    
    testing_results=pd.DataFrame(prepped_trip_frame[:len(trip_results)])
    testing_results['results'] = pd.Series(trip_results)
    testing_results['mean mi/kwh'] = testing_results['results'].apply(lambda x: pd.Series(x).mean())
    testing_results['std mi/kwh'] = testing_results['results'].apply(lambda x: pd.Series(x).std())
    testing_results['med mi/kwh'] = testing_results['results'].apply(lambda x: pd.Series(x).median())
    return testing_results
    