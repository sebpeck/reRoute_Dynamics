'''
Instance_Tools.py
S. Peck

Instance_Tools.py contains methods used in the processing and handling of individual points on a route and other
miscellaneous methods in generating a trip.

Methods:
generate_riders() - method used to generate a series of ridership changes at each stop based on the
expected ridership of that series.
check_hit_signal() - method to return if a signal light has been hit or not based on random chance.
determine_stop_type() - method used to determine what kind of stop a position is based on ridership, signals, and signs.
get_stop_type() - method to create an array of stop type booleans for a given position.
get_distances_to_stop() - method to determine how much distance is between the current index and the next stop.

'''

import numpy as np
import random
import pandas as pd


def generate_riders(n_stops, mean_ridership, seed=None):
    '''
    generate_riders() uses a number of stops, and mean ridership,
    and then generates a list of riderhsip changes by stops, such that
    the ridership is as close to the ceiling of the mean as possible.
    
    Params:
    n_stops - int of number of stops.
    mean_ridership - float of the mean ridership for the trip
    seed - int, default None, used for the seed of the randomness.
    
    Returns:
    list of randomly generated ridership changes at each stop.
    '''
    
    # Check if the seed isn't None/
    if seed is not None:
        
        # Implement the seed.
        random.seed(seed)
    
    if n_stops > 1:
        # create a list of indexes
        rider_list_indexes = list(np.arange(0, n_stops, 1))

        # create a list of zeros based on number of stops
        riders_on = [0]*n_stops

        # for each rider in the mean, randomly select a stop they get on at.
        for i in range(int(np.ceil(mean_ridership))):
            riders_on[random.choice(rider_list_indexes[:-1])] += 1
            
        # tally for current riders.
        current_riders = 0

        # list to store ridership changes
        rider_changes = []

        # loop through each stop in the ridership on changes
        for stop in riders_on:

            # a number between 0 and the current number of riders disembark.
            riders_off = random.randrange(0, current_riders+1)

            # The current riders is updated for the net change in ridership.
            current_riders += stop - riders_off

            # the ridership change is calculated.
            rider_change = stop - riders_off

            # the change is appended.
            rider_changes.append(rider_change)


        # if not everyone is off by the last stop, shove em off.
        if current_riders != 0:
            rider_changes[-1] -= current_riders

        # return the rider changes.
        return rider_changes
    
    else:
        return [0]


def check_hit_signal(stoplight_chance=.541666, seed=None):
    '''
    check_hit_signal() takes a chance to hit a yellow/red,
    and randomly generates a flag if the light is hit (a stop) or not.
    
    Params:
    stoplight_chance - float, the fraction of time a stoplight will be hit. defualt based on 120s cycle, 55s g, 65 red/yellow. https://wsdot.wa.gov/travel/operations-services/traffic-signals
    seed - int, default none, used for the seed of randomness.
    
    Returns:
    an int, 0 or 1, depending if the light is green or red.
    '''
    
    # Check if the seed isn't None/
    if seed is not None:
        
        # Implement the seed.
        random.seed(seed)
        
    # generate a number between 1 and 100, if it's lower than the chance, it's a stop.
    hit = int(random.randrange(0, 101)<stoplight_chance*100)
    
    # return the generated hits.
    return hit


def determine_stop_type(rider_changes, hit_signals, signs):
    '''
    determine_stop_type() is used to determine the stop types based on
    ridership, signals, and signs, and ends.
    
    Params:
    rider_changes - iterable of ridership changes for a route
    hit_signals - iterable of signals for a route
    signs - iterable of signs for a route.
    
    Returns:
    a list of lists containing 0 or 1 for each stop type, with index of the sub-lists
    corresponding to the passed stop types.
    '''
    
    
    stop_data =  pd.concat([
                 pd.Series(rider_changes),
                 pd.Series(hit_signals),
                 pd.Series(signs)], axis=1)
    stopinfo = list(stop_data.apply(lambda x: get_stop_type(x[0], x[1], x[2], 0), axis=1))
    
    # make sure the ends have stops.
    stopinfo[0][3] = 1
    stopinfo[-1][3] = 1
    
    return stopinfo
    

def get_stop_type(ridership, signal, sign, end):
    '''
    get_stop_type() takes a ridership, signal, and sign point,
    and converts to a list of 0's and 1's depending on if each is
    valid or not.
    '''
    bus_stop_flag = False
    signal_flag = False
    sign_flag = False
    end_flag = False
    
    if ridership != 0:
        bus_stop_flag = True
    if signal != 0:
        signal_flag = True
    if sign != 0:
        sign_flag = True
    if end != 0:
        end_flag = True
    
    return [int(bus_stop_flag), int(signal_flag), int(sign_flag), int(end_flag)]


def get_distances_to_stops(stop_types, traveled_distance):
    '''
    get_distances_to_stops() takes an iterable of stop types, and the corresponding cumulative travel distances,
    and then generates an iterable of distance to the next stop of any type.
    
    Params:
    stop_types - iterable of stop types as generated by determine_stop_type()
    traveled_distance - iterable of cumulative traveled distance.
    
    Returns:
    iterable of the distance to the next stop based on current locaton.
    '''
    
    # join the two iterables
    df = pd.concat([pd.Series(stop_types), pd.Series(traveled_distance)], axis=1)
    
    # check if the stop types indicate any stop
    df[0]= df.apply(lambda x: sum(x[0]), axis=1) > 0
    
    # filter to only stops
    stops_and_distances = df[df[0] == True]
    
    # first position is a stop and should be marked as such
    dx_to_stop = [0]

    # loop through the remaining
    for i in range(1, len(stops_and_distances)):
        
        # get the current and last index of the filtered data
        current_index = stops_and_distances.index[i]
        last_index = stops_and_distances.index[i-1]
        
        # Get the data range between the indexes
        current_range = df.iloc[last_index+1:current_index+1][1]
        
        # Get the stop distance, (the last of the current range
        stop_distance = current_range.iloc[-1]
        
        # get the difference between the stop distance and each point
        current_range = abs(stop_distance-current_range)
        
        # extend the distance list with the recalculation
        dx_to_stop.extend(list(current_range))
        
    # return the data
    return dx_to_stop