'''
Trip_Simulator.py
'''
import Instance_Tools as it
import pandas as pd
import Object_Params as op
import Physics_Engine as pe
import Geography_Tools as gt
from alive_progress import alive_bar

def simulate_trip(route, trip=op.Trip(), bus=op.Bus(), ESS=op.ESS()):
    
    # make copies so there isnt accidental overwriting of params.
    trip = trip.copy()
    bus = bus.copy()
    ESS = ESS.copy()
    
    # take the stops and convert such that -1 is 0.
    rider_series = pd.Series(route.stops)
    rider_series = rider_series + 1

    # get total number of stops by filtering out the 0's
    n_stops = len(rider_series[rider_series != 0])

    # set the positions with stops to ridership changes instead using generate_riders
    rider_series[rider_series != 0] = it.generate_riders(n_stops, mean_ridership = trip.m_riders, seed=trip.seed)

    # convert the ridership data to a list.
    rider_changes = list(rider_series)

    # take the signal data and convert invalid points to 0
    signal_series = pd.Series(route.signals)
    signal_series[signal_series == -1] = 0

    # multiply the signal values by a boolean of if the stoplight has been hit or not.
    signal_series.apply(lambda x: x*it.check_hit_signal(stoplight_chance = trip.chance_sig, seed=trip.seed))


    # convert the data to a list.
    signals_hit = list(signal_series)

    # get the stop type arrays for each point
    stop_types = it.determine_stop_type(rider_changes, signals_hit, pd.Series(route.signs)+1)

    # use the stop type arrays and running distance to get the distance from each point to its next stop.
    stop_distances = it.get_distances_to_stops(stop_types, route.cum_d_X)

    # Convert stop distances to meters
    stop_distances = list(pd.Series(stop_distances)*1000)

    # Get the grades
    grades = route.grades

    # get the grade forces
    grade_forces = list(pd.Series(grades).apply(lambda x: pe.calculate_grade_force(x,bus.mass, bus.Cf)))

    # Get the changes in travel distance, in meters
    dxs = list(pd.Series(route.d_X)*1000)

    # get the speed limits, in m/s, and affect it by the trip's traffic parameter
    limits = list(pd.Series(route.limits)*1000*((1-(trip.traffic/2)))*(1+bus.f_a))

    # get the bearing angle from headings/
    bearings = list(pd.Series(route.bearings).apply(lambda x: gt.heading_to_angle(x)))

    # get the geodesic distances in meters
    geodes_dists =list(pd.Series(route.dx)*1000)

    # get the elevations, in km
    elevations = route.elevation
    
    geometry = route.geometry
    
    wind_angle = gt.heading_to_angle(trip.wind_bearing)

    # create an empty dict to store the result structure
    initialize_result = {'type':'initialize',
                       'v_f':0,
                       'dt':0,
                       'P':0,
                       'dx':0,
                       'grade':0,
                       'limit':0,
                       'stop_clf':[True, True, True, True],
                        'gdx':0,
                        'elevation':0,
                        'dx_to_next':0,
                        'b_dx':0,
                        'geometry':0}

    # Initialize the running data with the initialize dict.
    running_data = [initialize_result]

    # Loop through each point. 
    with alive_bar(len(stop_distances)-1) as bar:
        for i in range(len(stop_distances)-1):
            #print(running_data[-1])
            v = running_data[-1]['v_f']

            # Query point data
            dx_to_next_stop = stop_distances[i]
            grade = grades[i]
            grade_force = grade_forces[i]
            dx = dxs[i]
            limit = limits[i]
            bearing = bearings[i]
            stop_classes = stop_types[i+1]

            rider_change = rider_changes[i]
            geodes_dist = geodes_dists[i]
            elevation = elevations[i]
            position = geometry[i]

            # Calculate the wind force
            wind_force = pe.calculate_wind_force(bearing,
                                                 v,
                                                 wind_angle,
                                                 trip.v_wind,
                                                 trip.p_air,
                                                 bus.Cd,
                                                 bus.area)

            # Calculate the braking distance
            braking_distance_data = pe.get_braking_distance(v,
                                                            bus.mass,
                                                            grade_force,
                                                            wind_force,
                                                            bus.a_br,
                                                            bus.f_br,
                                                            bus.f_i,
                                                            bus.dmax)

            # Parse braking information
            braking_distance = braking_distance_data['dx']
            adjusted_braking_factor = braking_distance_data['bf']

            # Conditional booleans
            above_limit = (v > (limit + trip.MOE*limit))
            below_limit = (v < (limit - trip.MOE*limit))
            stopped = (v < trip.stop_margin)
            stop_upcoming = (braking_distance+limit*bus.dt_max >= dx_to_next_stop )

            # check if the distance to the next stop is not 0
            if dx_to_next_stop != 0:
                braking_distance_data = pe.get_braking_distance(v,
                                                            bus.mass,
                                                            grade_force,
                                                            wind_force,
                                                            bus.a_br,
                                                            bus.f_br,
                                                            bus.f_i,
                                                            dx_to_next_stop)
                    # Parse braking information
                braking_distance = braking_distance_data['dx']
                adjusted_braking_factor = braking_distance_data['bf']


            # Set up a dict for the true result.
            true_result = {'type':'NULL',
                           'v_f':None,
                           'dt':None,
                           'P':None,
                           'dx':dx,
                           'grade':grade,
                           'limit':limit,
                           'stop_clf':stop_classes,
                           'gdx':geodes_dist,
                           'elevation':elevation,
                           'dx_to_next':dx_to_next_stop,
                           'b_dx':(braking_distance, adjusted_braking_factor),
                           'geometry':position}

            # Check if there's a stop upcoming and if the bus is still moving
            if stop_upcoming and not stopped:
                #brake to stop.
                result = pe.brake(v,
                                  bus.mass,
                                  dx,
                                  grade_force,
                                  wind_force,
                                  bus.a_br,
                                  adjusted_braking_factor,
                                  bus.f_i,
                                  bus.dmax)



                # insert the results.
                true_result['type'] = 'stp_brk'
                true_result['v_f'] = result['v_f']
                true_result['dt'] = result['dt']
                true_result['P'] = result['P']
                running_data.append(true_result.copy()) #<-- This has to use the copy, otherwise it will change prev. values

                # if the speed is within the stop margin, stop the bus.
                if (result['v_f'] < trip.stop_margin):
                    result['v_f'] = 0
                    # Calculate the stop time based on the class.
                    stop_time = stop_classes[0]*trip.t_stop + stop_classes[1]*trip.t_sig + stop_classes[2]*trip.t_sign + stop_classes[3]*trip.t_end  # ridership, signal, sign, end

                    # Adjust the mass of the bus
                    bus.mass += rider_change*trip.m_pass

                    # Add the results.
                    true_result['type'] = 'rest'
                    true_result['v_f'] = 0
                    true_result['dt'] = stop_time
                    true_result['P'] = 0
                    true_result['dx'] = 0
                    tmp_storage = true_result['stop_clf']

                    running_data.append(true_result.copy())

            # otherwise, if the bus is not coming up to a stop:
            else:
                if stopped:
                    #accelerate

                    # Get a list of the results for accelerate (done like this as accelerate *may*
                    # give a number of dicts depending on how development goes
                    results = [pe.accelerate(v,
                                           bus.mass,
                                           dx,
                                           grade_force,
                                           wind_force,
                                           bus.a_prof,
                                           bus.a_br,
                                           bus.f_br,
                                           bus.f_i,
                                           bus.P_max,
                                           bus.a_max,
                                           bus.dt_max)]

                    # Append the results
                    for result in results:
                        true_result['type'] = 'ac_from_0'
                        true_result['v_f'] = result['v_f']
                        true_result['dt'] = result['dt']
                        true_result['P'] = result['P']
                        #print('ac_from_0:',result)
                        running_data.append(true_result.copy())
                    
                    #print("v, dx, grade_force, wind_force = {}, {}, {}, {}".format(v, dx, grade_force, wind_force))

                elif below_limit and not stopped:
                    #accelerate
                    results = [pe.accelerate(v,
                               bus.mass,
                               dx,
                               grade_force,
                               wind_force,
                               bus.a_prof,
                               bus.a_br,
                               bus.f_br,
                               bus.f_i,
                               bus.P_max,
                               bus.a_max,
                               bus.dt_max)]
                    for result in results:
                        true_result['type'] = 'ac_below'
                        true_result['v_f'] = result['v_f']
                        true_result['dt'] = result['dt']
                        true_result['P'] = result['P']
                        #print('ac_below:',result)
                        running_data.append(true_result.copy())

                elif above_limit and not stopped:
                    #brake

                    result = pe.brake(v,
                                  bus.mass,
                                  dx,
                                  grade_force,
                                  wind_force,
                                  bus.a_br,
                                  adjusted_braking_factor,
                                  bus.f_i,
                                  bus.dmax)
                    true_result['type'] = 'br_above'
                    true_result['v_f'] = result['v_f']
                    true_result['dt'] = result['dt']
                    true_result['P'] = result['P']
                    #print('br_above:',result)
                    running_data.append(true_result.copy())

                else:
                    #maintain
                    result = pe.maintain(v,
                                         bus.mass,
                                         dx,
                                         grade_force,
                                         wind_force,
                                         bus.a_br,
                                         bus.f_br,
                                         bus.f_i, 
                                         bus.P_max)
                    true_result['type'] = 'main'
                    true_result['v_f'] = result['v_f']
                    true_result['dt'] = result['dt']
                    true_result['P'] = result['P']
                    #print('main:',result)
                    running_data.append(true_result.copy())

                #print(running_data[-1])
            bar()

    # Not converted to ESS yet
    return running_data.copy()


if __name__ == '__main__':
    
    # load the arguments
    parser = argparse.ArgumentParser(description='Simulate a bus trip.')
    parser.add_argument('-r',"--route", help="Path to a Route export json from Geography_Tools.")
    parser.add_argument('-t',"--trip", help="Path to a trip parameter export txt from Object_Params.")
    parser.add_argument('-b',"--bus", help="Path to a bus parameter export txt from Object_Params.")
    parser.add_argument('-e',"--ESS", help="Path to a ESS parameter export txt from Object_Params.")
    parser.add_argument('-o', "-output", help="Output filepath.")
    parser.add_argument('-n', "-name", help='Output filename.')
    
    # parse the arguments
    args = parser.parse_args()
    
    #  Check for the optional params
    if type(args.trip) == type(None): trip_obj = OP.Trip()
    else: trip_obj = op.load_trip_params(args.trip)
    
    if type(args.bus) == type(None): bus_obj = OP.Bus()
    else: bus_obj = op.load_bus_params(args.bus)
    
    if type(args.ESS) == type(None): ESS_obj = OP.ESS()
    else: ESS_obj = op.load_ESS_params(args.ESS)

    # get the path params
    output_path = args.output
    output_name = args.name
    
    # Check the required params
    if type(args.route) == type(None): route_obj = gt.load_from_json(input("[ Trip_Simulator.py ] Route filepath: "))
    if type(output_path) == type(None): output_path = input("[ Trip_Simulator.py ] Output path: ")
    if type(output_name) == type(None): output_name = input("[ Trip_Simulator.py ] Output name: ")

    # run the simulation
    results = simulate_trip(route_obj, trip_obj, bus_obj, ESS_obj)
    
    # Save the data as a feather.
    pd.DataFrame(results).to_feather('{}/{}.feather'.format(output_path, output_name))
    
    
    