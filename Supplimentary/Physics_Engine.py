'''
Physics_Engine.py
'''

import numpy as np
import pandas as pd
import Geography_Tools as gt

DEFAULT_AIR_DENSITY = 1.2 #kg/m^3, approximate of 20 deg C.
DEFAULT_WIND_SPEED = 4 #m/s Typical for seattle
DEFAULT_WIND_BEARING = gt.heading_to_angle('SE') # Typical for seattle
DEFAULT_DRAG_COEFF = 0.6 # From Aggregate, Abdelaty & Mohamed []
DEFAULT_FRONTAL_AREA = 2.6 * 3.3 # meters, width times height New Flyer Xcelsior Charge []

DEFAULT_COEFF_F = .01 # Abdelaty & Mohamed

DEFAULT_BR_ACCEL = 1.5 #m/s^2 https://www.apta.com/wp-content/uploads/APTA-BTS-BC-RP-001-05_Rev1.pdf <-- Possible source, handbrake road minimum is ~1.5
DEFAULT_BR_FACTOR = .5 # Half of braking accel.
DEFAULT_I_FACTOR = 1.1 # intertial factor, Asamer Et. Al []
DEFAULT_MAX_DIST = 304.8 # meters, approx based on offramp length for I-5

DEFAULT_MAX_POWER = 160000 # Watts

DEFAULT_A_PROF = pd.read_csv("./Data/Accel_Prof/acceleration.csv", header=None) # Default acceleration profile.
DEFAULT_A_PROF[1] = DEFAULT_A_PROF[1].apply(lambda x: x*9.81)
DEFAULT_MAX_ACC = .4 #m/s^2 the asymptotic acceleration a profile will achieve.
DEFAULT_MAX_DT = 1 #second, the timestep in the profile corresponding to the max acceleration.



def calculate_wind_force(bearing,
                         speed,
                         wind_bearing = DEFAULT_WIND_BEARING,
                         wind_speed = DEFAULT_WIND_SPEED,
                         air_density = DEFAULT_AIR_DENSITY,
                         drag_coeff = DEFAULT_DRAG_COEFF,
                         frontal_area = DEFAULT_FRONTAL_AREA):
    '''
    calculate_wind_force() is used to calculate the force exerted against
    the bus's direction of travel by the wind.
    
    Params:
    bearing - compass heading, in degrees (float), the bus is travelling in.
    speed - current bus velocity. (float)
    wind_bearing - compass heading, in degrees (float), the wind is in. Default of SE directon
    wind_speed - current wind velocity in its direction of travel. (float) Default of 4 m/s
    air_density - current air density. (float), default 1.2 kg/m^3
    drag_coeff - drag coefficient of the bus. (float), default .6.
    frontal_area - frontal area of bus.  (float), default 2.6*3.3 m^2
    
    Returns:
    float representation of force exerted against bus's direction
    of travel. Negative means the bus is being accelerated by the wind.
    '''
    
    # convert bearing to angle between -180 and 180:
    centered_bearing = bearing - 180
    centered_wind_bearing = wind_bearing - 180
    
    # get wind angle relative to bus:
    relative_wind_angle = centered_wind_bearing - centered_bearing
    
    # get velocity of wind relative to bus:
    relative_wind_v = -wind_speed*np.cos(np.radians(relative_wind_angle))
    
    air_force = sign(relative_wind_v)*drag_coeff*frontal_area*air_density/2*(speed - relative_wind_v)**2
    
    return air_force


def sign(num):
    '''
    sign checks the sign of a number.
    '''
    return -1 if num < 0 else 1

    
def calculate_grade_force(grade, mass, f_coeff=DEFAULT_COEFF_F):
    '''
    calculate_grade_force takes a grade of a slope, a mass, and friction
    coefficient, and calculates the force due to gravity against the direction of
    travel of the mass. Negative value means the object is being accelerated, 
    rather than decellerated.
    
    Params:
    grade - grade angle of a point, as a percent float.
    mass - mass of object as a float, in kg
    f_coeff - the coefficient of friction. default .01.
    
    Returns:
    net force due to gravity as float, in Newtons
    '''
    
    # get the angle of the grade
    grade_angle = np.arctan(grade/100)
    
    # calculate the hill force
    hill_force = np.sin(grade_angle) * 9.81 * mass
    
    # calculate the frictional force (always positive and opposing direction of travel)
    fric_force = np.cos(grade_angle) * 9.81 * mass * f_coeff
    
    # return the net force on the object at that point.
    return hill_force + fric_force


def get_braking_distance(velocity,
                         mass,
                         grade_force,
                         wind_force,
                         braking_acceleration = DEFAULT_BR_ACCEL,
                         braking_factor = DEFAULT_BR_FACTOR,
                         intertial_factor = DEFAULT_I_FACTOR,
                         max_distance = DEFAULT_MAX_DIST
                         ):
    
    '''
    get_braking_distance() is used to determine how far an object needs to stop.
    
    Params:
    velocity - object's initial velocity as a float, in m/s.
    mass - mass of the object in kg.
    grade_force - force, in N, due to the grade of the slope the object experiences.
    wind_force - force, in N, experienced by the object due to wind resistance.
    braking_acceleration - the maximum acceleration due to braking, as float. Default 1.5 m/s^2
    braking_factor - float between 0 and 1 representing how much of the max braking acceleration
                     is used. Default .5.
    inertial_factor - float representing how inertia affects the acceleration of the bus. Default 1.1
    max_distance - float representing the ideal maximum distance the bust takes to stop. Default 304.8 m.
    
    Returns:
    a dict containing braking distance dx, braking factor bf, and decelleration rate ad.
    
    '''
    # in other words, the force being applied against the bus's direction of travel before braking
    net_external_force = grade_force + wind_force # positive means decellerating the bus.
    
    # calculate an initial pass of decelleration
    rate_decell = (braking_acceleration*braking_factor*intertial_factor) + net_external_force/mass

    # check if the current deceleration rate is higher than the 'maximum that's expected'
    while (rate_decell > braking_acceleration*intertial_factor) and braking_factor > 0.01:

        # reduce braking factor
        braking_factor -= .0001
        
        # recalculate the decelleration
        rate_decell = (braking_acceleration*braking_factor*intertial_factor) + net_external_force/mass
    
    
    # Check to make sure the object is actually decelerating to be able to stop over the max distance.
    while (braking_factor<1.001) and rate_decell < velocity**2/(2*max_distance):
        
        # if it isn't, up the braking factor.
        braking_factor += .0001
        
        # update the rate of deceleration.
        rate_decell = (braking_acceleration*braking_factor*intertial_factor) + net_external_force/mass

            
    # calculate the braking distance.
    braking_dist = velocity**2 / (2*rate_decell)
    
    # return the braking distance.
    return {'dx':braking_dist, 'bf':braking_factor, 'ad':rate_decell}


def brake(velocity,
          mass,
          travel_distance,
          grade_force,
          wind_force,
          braking_acceleration = DEFAULT_BR_ACCEL,
          braking_factor = DEFAULT_BR_FACTOR,
          inertial_factor = DEFAULT_I_FACTOR,
          max_distance = DEFAULT_MAX_DIST):
    '''
    brake() is used to determine the velocity change, time change, and power consumption of braking
    a mass at a given velocity over a set distance.
    
    Params:
    velocity - object's initial velocity as a float, in m/s.
    mass - mass of the object in kg.
    travel_distance - distance to brake over, in meters.
    grade_force - force, in N, due to the grade of the slope the object experiences.
    wind_force - force, in N, experienced by the object due to wind resistance.
    braking_acceleration - the maximum acceleration due to braking, as float. Default 1.5 m/s^2
    braking_factor - float between 0 and 1 representing how much of the max braking acceleration
                     is used. Default .5.
    inertial_factor - float representing how inertia affects the acceleration of the bus. Default 1.1
    max_distance - float representing the ideal maximum distance the bust takes to stop. Default 304.8 m.
    
    Returns:
    a dict containing final velocity v_f, time change dt, power P, and braking factor bf.
    '''
    vf = 0
    t = 99999
    P = 99999
    bf = braking_factor
    
    if velocity != 0:

        # initial rate of deceleration
        rate_decell = (braking_acceleration * inertial_factor * braking_factor) + (grade_force + wind_force)/mass
        # first pass for vf
        pre_root = -2*rate_decell*travel_distance + velocity**2

        if pre_root >0:
            vf = np.sqrt(-rate_decell*2*travel_distance+(velocity**2))
        else:
            vf = 0

        # calculate time
        t = (velocity-vf)/rate_decell

        # use the bus's acceleration without externalities to calculate power needed
        P = mass*-(braking_acceleration * inertial_factor * braking_factor)*travel_distance/t

        # Check if velocity is 0 to re-do the power
        if vf == 0:
            P = 0
    else:
        vf = 0
        t = 1
        P = 0

        # create a dict of results
    results = {'v_f':vf,
               'dt':t,
               'P':P,
               'bf':braking_factor}
    
    # return the results
    return results
        

def maintain(velocity,
             mass,
             travel_distance,
             grade_force,
             wind_force,
             braking_acceleration = DEFAULT_BR_ACCEL,
             braking_factor = DEFAULT_BR_FACTOR,
             inertial_factor = DEFAULT_I_FACTOR,
             max_power = DEFAULT_MAX_POWER):
    '''
    maintain() is used to determine the velocity change, time change, and power consumption of maintaining
    a mass at a given velocity over a set distance.
    
    Params:
    velocity - object's initial velocity as a float, in m/s.
    mass - mass of the object in kg.
    travel_distance - distance to brake over, in meters.
    grade_force - force, in N, due to the grade of the slope the object experiences.
    wind_force - force, in N, experienced by the object due to wind resistance.
    braking_acceleration - the maximum acceleration due to braking, as float. Default 1.5 m/s^2
    braking_factor - float between 0 and 1 representing how much of the max braking acceleration
                     is used. Default .5.
    inertial_factor - float representing how inertia affects the acceleration of the bus. Default 1.1
    max_power - float representing the maximum motor power. Default 160000 W.
    
    Returns:
    a dict containing final velocity v_f, time change dt, and power P.
    
    '''
    
    # Calculate the counter-acceleration to counteract all external forces
    a_counter = (grade_force + wind_force)/mass
    
    # set variables
    P_max = max_power
    
    # Current power is the power needed to counteract all external forces
    P_current = mass*a_counter*velocity
    
    # P_brake_max is the maximum power the bus can output from braking
    P_brake_max = -mass*velocity*(1*braking_acceleration*inertial_factor)
    
    # others
    v_f = None
    v_i = velocity
    dt = None
    P = None
    m = mass
    dx = travel_distance
    
    # Check if there is enough motor power
    #print(P_max, P_current, P_brake_max)
    if P_max >= P_current > 0:
        #print("Enough motor power")
        v_f = v_i
        dt = v_f*dx
        P = P_current
    
    # Check if there isn't enough motor power to overcome.
    elif P_current > P_max:
        #print("Not enough motor power")
        da = (P_current - P_max)/(m*v_i)
        v_f = np.sqrt(v_i**2 - 2*dx*da)
        P = P_max
        dt = (v_i - v_f) / da
    
    # Check if there's enough braking power to overcome.
    elif 0 >= P_current > P_brake_max:
        #print("Enough braking power")
        v_f = v_i
        dt = v_f*dx
        P=P_current
    
    # Check if there's not enough braking power to overcome.
    elif P_brake_max >= P_current:
        #print("Not enough braking power")
        da = -abs(P_current-P_brake_max) / (m*v_i)
        v_f = np.sqrt(v_i**2 + 2*dx*da)
        P = P_brake_max
        dt = (v_i - v_f) / da
    
    # Otherwise, error time!
    else:
        print("You shouldn't be here")
        print('mass {}, a_counter {}, velocity {}, braking_Factor {}, braking_accel {}, i_factor {}'.format(mass, a_counter, velocity, braking_factor, braking_acceleration, inertial_factor))
        print("Pmax {}, pcurrent {}, pmin {}".format(P_max, P_current, P_brake_max))
        
    # turn the results into a dict
    results = {'v_f':v_f,
               'dt':dt,
               'P':P}
    
    # return the results
    return results


def accelerate(velocity,
               mass,
               travel_distance,
               grade_force,
               wind_force,
               raw_a_prof = DEFAULT_A_PROF,
               braking_acceleration = DEFAULT_BR_ACCEL,
               braking_factor = DEFAULT_BR_FACTOR,
               inertial_factor = DEFAULT_I_FACTOR,
               max_power = DEFAULT_MAX_POWER,
               max_acc = DEFAULT_MAX_ACC,
               max_timestep = DEFAULT_MAX_DT):
    '''
    accelerate() is used to determine the velocity change, time change, and power consumption of accelerating
    a mass at a given velocity over a set distance, with a given acceleration profile.
    
    Params:
    velocity - object's initial velocity as a float, in m/s.
    mass - mass of the object in kg.
    travel_distance - distance to brake over, in meters.
    grade_force - force, in N, due to the grade of the slope the object experiences.
    wind_force - force, in N, experienced by the object due to wind resistance.
    raw_a_prof - an acceleration profile, with the 0 index column being cumulative time,
                 and the 1 being acceleration at that time (in m/s^2)
    braking_acceleration - the maximum acceleration due to braking, as float. Default 1.5 m/s^2
    braking_factor - float between 0 and 1 representing how much of the max braking acceleration
                     is used. Default .5.
    inertial_factor - float representing how inertia affects the acceleration of the bus. Default 1.1
    max_power - float representing the maximum motor power. Default 160000 W.
    
    Returns:
    a dict containing final velocity v_f, time change dt, and power P.
    '''
    
    a_counter = (grade_force + wind_force)/mass
    #print(a_counter)
    
    # set variables to adjust later.
    P_max = max_power
    dt = 0
    v_f = 0
    P = 0
    a_max=0
    
    # Assume the bus can always start up. If the velocity is 0,
    if velocity ==0:
        
        # get the maximum acceleration as the first value + enough to counter.
        a_max = a_counter + raw_a_prof[1].iloc[0]
        
        # while the power is greater than the max power, reduce the max acceleration
        while a_max*mass*travel_distance/max_timestep > P_max:
            a_max -= .0001
    else:
        # Otherwise, set the max acceleration based on max power.
        a_max = P_max/(mass*velocity)
    
    # get a copy of the profile
    a_prof_ext = raw_a_prof.copy()
    
    # get the profile of the acceleration.
    a_prof = pd.concat([raw_a_prof]).reset_index(drop=True)
    
    # Set the last index of the acceleration profile as the max timestep and accel.
    for i in range(40):
        
        # Also, extend the profile with the timestep and acceleration.
        a_prof.at[-1, 0] = max_timestep + list(a_prof[0])[-1]
        a_prof.at[-1, 1] = max_acc
        a_prof = a_prof.reset_index(drop=True)
    
    #Convert the acceleration to needed acceleration.
    needed_a = a_prof[1].apply(lambda x: (x+a_counter)*inertial_factor)
    
    # convert times to time steps
    a_prof[0] = a_prof[0].diff().shift(-1)
    a_prof.at[len(a_prof)-1, 0] = max_timestep
    

    # get the regulator velocities
    reg_vel = a_prof[0]*a_prof[1].cumsum()
    
    # using the regulator velocities, get the starting index:
    starting_index=(reg_vel.apply(lambda x: abs(x-velocity)).sort_values().index)[0]
    
    # P = m*a_engine*(v+dt*a_engine)
    powers = mass*needed_a*(velocity + a_prof[0]*needed_a)
    powers = powers.clip(0, P_max)
    
    # Quadratic formula to solve for possible a from power equation
    possible_a = (-velocity+np.sqrt(velocity**2 + 4*a_prof[0]*powers/mass))/(2*a_prof[0])
    possible_a = possible_a.fillna(possible_a.min()) # USE THIS FOR CALCULATING POWER
    
    # a_net = a_engine/f_i - a_counter
    true_a = possible_a/inertial_factor - a_counter # USE THIS FOR CALCULATING TIME AND DISTANCE AND VELOCITY
    clipped_true_a = true_a.clip(0)

    # dv = a * dt
    vels = clipped_true_a * a_prof[0]
    
    # Create a profile dataframe using the existing data so as to index through it. 
    a_prof = pd.concat([a_prof[0], possible_a, clipped_true_a, powers, vels, (vels * a_prof[0]).cumsum()], axis=1)
    a_prof.columns = ['dt','int_a','net_a','power', 'dv', 'dx']
    a_prof['cum_dv'] = a_prof['dv'].cumsum()
    a_prof['cum_dx'] = a_prof['dx'].cumsum()
    a_prof['cum_dx'] = a_prof['cum_dx'] - a_prof['cum_dx'].iloc[starting_index] - travel_distance
    
    # Get the closest index to the end of the acceleration profile.
    closest_end_index = list(abs(a_prof['cum_dx']).sort_values().index)[0]
    
    #Check if the starting and ending indexes are the same 
    if (starting_index == closest_end_index):
        #print('same')
        # get the data from the index
        ea_prof = a_prof.iloc[starting_index]
        
        # Calculate final velocity
        v_f = np.sqrt(2*travel_distance*ea_prof['net_a'] + velocity**2)
        
        # Calculate change in time
        dt = -(velocity - v_f) / ea_prof['net_a']
        
        # Calculate power.
        P = mass*ea_prof['int_a']*travel_distance/dt

    # if the forces are too hard to overcome, the starting index will
    # be greater than the end index
    elif(starting_index > closest_end_index):
        #print('overforce')
        # get the acceleration from the true net accels.
        a=true_a.iloc[starting_index]
        # calculate the velocity
        v_f = np.sqrt(2*travel_distance*a + velocity**2)
        #print(travel_distance, a, velocity)
        
        # Calculate time change.
        dt = abs((velocity-v_f)/a)
        
        # use the max power.
        P = P_max
    
    # otherwise, 
    else:
        #print('range')
        #print("multi")
        # get the selection range of data
        ea_prof = a_prof.iloc[starting_index:closest_end_index]
        
        # initialize the final velocity
        v_f = velocity
        
        # intialize the total energy
        en=0
        dx = 0
        # iterate through the range and adjust the
        # velocity, time, and energy accordingly
        for col, row in ea_prof.iterrows():
            v_f = v_f+row['dv']
            i_dt = row['dt']
            dt += i_dt
            en += row['power']*i_dt
            dx += .5*(v_f + v_f-row['dv'])*dt
            #print(dx)
            
        #print(dx)
        # calculate the power.
        P = en/dt
        
    
    # Get the results.
    results = {'v_f':v_f,
               'dt':dt,
               'P':P}
    # return the results
    return results