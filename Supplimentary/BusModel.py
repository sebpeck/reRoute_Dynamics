import pandas as pd
import numpy as np

class BusModel:
    def __init__(self,
                 acceleration_profile, # two column dataframe containing time(s), accel(m/s^2)
                 raw_mass = 13041, # kilograms, New Flyer Xcelsior Charge []
                 bus_width = 2.6, # meters, New Flyer Xcelsior Charge []
                 bus_height = 3.3, # meters, New Flyer Xcelsior Charge []
                 drag_coeff = 0.6, # From Aggregate, Abdelaty & Mohamed []
                 factor = 1.1, # intertial factor, Asamer Et. Al []
                 fric_coeff = .01, # Abdelaty & Mohamed []
                 motor_eff = .916, # unknown source, not listed in paper -- tourque and speed depemdamt, next level would be based on tourqe
                 invert_eff = .971, # unknown source, not listed in paper
                 max_power = 160, # kW, output New Flyer Xcelsior Charge []
                 regen = .6, # Gallet et al []
                 eff_aux = .89, # unknown source, Auxulliary system efficiency?
                 aux_load = 7000, # W, Abdelaty & Mohamed []
                 a_braking = -1.5, #m/s^2 # https://www.apta.com/wp-content/uploads/APTA-BTS-BC-RP-001-05_Rev1.pdf <-- Possible source, handbrake road minimum is ~1.5
                 final_a = .4, # m/s^2, defualt acceleration after profile finishes
                 max_velocity = 26.8224, # m/s, = 60 mph
                 maintain_acceleration = False, # boolean for if bus should maintain the last a in 
                                                # profile or use final_a for extrapolating new vals
                 num_starting_passengers = 0,
                 pass_ave_mass = 80, #kg, Gallet et Al. []
                
                 ):
        # raw bus characteristics -- Not all of these are used, but may be useful for the future
        self._empty_mass = raw_mass
        self._bus_width = bus_width
        self._bus_height = bus_height
        self._bus_front_area = self._bus_width * self._bus_height
        self._drag_coeff = drag_coeff
        self._i_factor = factor
        self._fric_coeff = fric_coeff
        self._motor_eff = motor_eff
        self._invert_eff = invert_eff
        self._max_power = max_power
        self._regen = regen
        self._eff_aux = eff_aux
        self._aux_load = aux_load
        
        # driving model characteristics
        self._a_braking = a_braking
        self._final_a = final_a
        self._final_a_boolean = maintain_acceleration
        self._max_v = max_velocity
        self._accel_profile_df = self.interpret_accel_prof(acceleration_profile)
        
        # variables based on position in route
        self._passengers = num_starting_passengers
        self._mass_p_pass = pass_ave_mass
        self._current_velocity = 0
        self._current_accel = 0
        self._distance_traveled = 0
        self._bus_status = "Stopped"
        self._current_mass = self._empty_mass + self._passengers * self._mass_p_pass
        
        
    def interpret_accel_prof(self, accel_prof):
        """
        interpret_accel_prof is used to format the acceleration
        profile so as to have velocity and distance traveled.
        
        Parameters:
        accel_prof: dataframe of time and acceleration values,
                    in seconds and m/s^2.
        
        Returns: 
        dataframe of time, acceleration, velocity, and distance
        in seconds and meters units.
        """
        # Generate a copy of the profile
        acc_prof = accel_prof.copy()
        
        # Get the changes in time
        acc_prof['dt'] = acc_prof['time[s]'].diff()
        
        # use mean_integral to get the velocity data
        acc_prof['vel.[m/s]'] = self.mean_integral(acc_prof['dt'], acc_prof['accel.[m/s^2]'])
        
        # use mean_integral to get the distance data
        acc_prof['dist.[m]'] = self.mean_integral(acc_prof['dt'], acc_prof['vel.[m/s]'])
        
        while (not (acc_prof['vel.[m/s]'].iloc[-1] >= self._max_v)):
            
            # Generate a new row dictionary from the columns
            new_row = dict.fromkeys(list(acc_prof.columns))
            
            # set the step time to be 1 second
            d_t = 1
            
            # calculate the next time point using the previous and dt
            next_time = acc_prof['time[s]'].iloc[-1] + d_t 
            
            # default to using the _final_a acceleration
            next_accel = self._final_a
            
            # if the parameters say to maintain acceleration,
            if self._final_a_boolean:
                # Keep the last acceleration
                next_accel = acc_prof['accel.[m/s^2]'].iloc[-1]
            
            # calculate the next velocity
            next_v = acc_prof['vel.[m/s]'].iloc[-1] + d_t*next_accel
            
            # calculate the next distance
            next_d = acc_prof['dist.[m]'].iloc[-1] + d_t * (next_v + acc_prof['vel.[m/s]'].iloc[-1]) / 2
            
            # generate the list of data
            data = [next_time, next_accel, d_t, next_v, next_d]
            
            # put the data into the dictionary
            for i in range(len(data)): new_row[list(new_row.keys())[i]] = data[i]
            
            # Append the new row on to the profile dataframe
            acc_prof.loc[len(acc_prof)] = new_row
            
        # drop the dt column
        #acc_prof = acc_prof.drop(columns = ['dt'])
        
        # return the acceleration profile dataframe
        return acc_prof
    
    
    def mean_integral(self, time_changes, value_series):
        """
        mean_integral takes in a series of dt (s) values, and
        a corresponding series of other values to be integrated
        with respect to time, and provides the result as a
        cumulative summing series. 

        Parameters:
        time_changes: a pandas series of dT values in units of seconds
        value_series: a pandas series of d[value] that has relation to time.

        Returns: 
        a pandas series of the integrated values.
        """
        int_prof = pd.DataFrame()
        int_prof['mean_val_btw'] = value_series.rolling(window=2).mean()
        int_prof['d[value]'] = int_prof['mean_val_btw'] * time_changes
        int_prof['integrated'] = int_prof['d[value]'].cumsum()

        return int_prof['integrated'].fillna(0)
    
    
    def get_braking_distance(self, velocity, braking_factor, ext_acc):
        '''
        get_braking_distance() takes in a velocity and the current external accelerations,
        and provides a calculation as to how far the bus would have to brake to reach a velocity of zero.
        
        Parameters:
        velocity: velocity, in m/s, the bus is travelling at
        ext_acc: the external accelerative force in m/s^2 the bus is experiencing.
        
        Returns:
        a distance, in meters, of how far the bus would have to brake to reach a velocity of zero.
        '''
        
        
        wind_accel = self.get_aerodynamic_drag(velocity, 0, 1.2)
        
        # Combine to get the net acceleration on the bus
        net_accel = -wind_accel - ext_acc
        a_bus = net_accel + self._a_braking*braking_factor*self._i_factor
        # Using the kinematics equation vf^2 = vi^2 + 2a(dX) to get
        # braking distance
        
        # ensure braking is not too harsh
        while (a_bus < self._a_braking):
            
            #add to the braking factor
            braking_factor -= .01
            
            # recalc the braking acceleration and bus acceleration
            braking_acc = self._a_braking*braking_factor*self._i_factor
            a_bus = net_accel + braking_acc
        
        # make sure that the braking factor is enough to actually slow the bus
        while (a_bus > -0.001 ) and (braking_factor <= 1):
            
            # Add to the braking factor
            braking_factor += .01
            
            # recalc the braking acceleration and bus acceleration
            braking_acc = self._a_braking*braking_factor*self._i_factor
            a_bus = net_accel + braking_acc
        
        dist = round((-velocity**2 / (2*(a_bus))), 5)
        
        # Check if the braking distance is not possible given the parameters
        if (dist < 0):
            dist = 2105 # Braking distance for mass of 25000, grade of -17 (seattle max), speed of 30 m/s, full brake
        return dist

    
    def get_aerodynamic_drag(self, bus_speed, wind_speed, air_density):
        '''
        get_aerodynamic_drag takes in wind speed and air density, and returns
        the acceleration, in m/s^2, they provide.
        
        Parameters:
        wind_speed: velocity, in m/s, of the wind.
        air_density: density, in kg/m^3, of the air.
        
        Returns:
        aerodynamic acceleration, in m/s^2.
        '''
        # Using the drag coefficient, bus frontal area, and passed parameters, calculate the acceleration of air due to drag
        air_drag = (self._drag_coeff * self._bus_front_area * (air_density/2) * (bus_speed - wind_speed)**2)/self._current_mass
        
        # Return said value.
        return air_drag # acceleration of air drag
    
    
    def get_inertial_acceleration(self):
        '''
        get_inertial_acceleration calculates the current acceleration of intertia the bus is undergoing.
        
        Parameters:
        N/A
        
        Returns:
        inertial acceleration of the bus, in Newtons.
        
        '''
        
        # Use the intertial factor, current mass, and current acceleration,
        # to get inertial acceleration.
        inertia = (self._i_factor * self._current_mass * self._current_accel) / self._current_mass
        return inertia # acceleration of inertia
    
    
    def update_riders(self, val, cat='change'):
        '''
        update riders takes in a value of how many riders get on (positive),
        or disembark (negative), and affects the bus mass and passenger count of the bus.
        
        Parameters:
        val: change in passengers as a float.
        cat: string, either 'change' or 'set', depending on if the number of passengers
             is being changed, or set. Default is 'change'.
        '''
        # Update passenger number, sum should never be < 0
        if (cat=='change'):
            self._passengers = self._passengers + val
        elif (cat == 'set'):
            self._passengers = val
        
        # Update total mass
        self._current_mass = self._empty_mass + self._passengers * self._mass_p_pass
        return None
    
    
    def brake_v3(self, dist, braking_factor, ext_acc):
        '''
        brake_v3() takes in a distance and braking factor and external acceleration,
        and determines the power and elapsed time for the bus to
        experience those conditions.
        
        Parameters:
        dist: distance to brake, in meters, as a float
        braking_factor: float representing how aggresive braking will be (max:1, min:0)
        ext_acc: the external acceleration, as a float, the bus
                 is experiencing due to gravity, friction, and elevation change.
        
        Returns:
        a touple containing the power [W] used and the elapsed time [s]. 
        '''

        # Get the current velocity and braking acceleration, 
        # using the factor to help determine how intense the braking is.
        v0 = self._current_velocity
        braking_acc = self._a_braking*braking_factor*self._i_factor

        # determine the external accelerations from intrinsic bus factors
        route_accel = ext_acc
        wind_accel = self.get_aerodynamic_drag(v0, 0, 1.2)

        # Generate the net acceleration the bus is experiencing, a_bus
        net_accel = (- wind_accel - route_accel)
        a_bus = net_accel + braking_acc
        
                # adjust braking to make sure it is less that .15 g's of acceleration
        while (a_bus < self._a_braking):
            braking_factor -= .01
            # recalc the braking acceleration and bus acceleration
            braking_acc = self._a_braking*braking_factor*self._i_factor
            a_bus = net_accel + braking_acc
        
        # make sure that the braking factor is enough to actually slow the bus
        while (a_bus > -0.001 ) and (braking_factor <= 1):
            
            # Add to the braking factor
            braking_factor += .01
            
            # recalc the braking acceleration and bus acceleration
            braking_acc = self._a_braking*braking_factor*self._i_factor
            
            a_bus = net_accel + braking_acc

            
        vf = 0
        # use kinematic equation to calculate final velocity,
        # Check to ensure calc is possible.
        if (v0**2 > np.abs(2*a_bus*dist)):
            vf = np.sqrt(v0**2 + 2*a_bus*dist)
            # if final velocity is slower than .1 m/s, set it to be zero
            if (vf < .1):
                vf = 0
        else:
            vf = 0


        #get the change in time from kinematics equations
        dt = ((vf - v0)/a_bus)

        # Calculate the power of the braking
        power_calc = braking_acc*self._current_mass*dist/dt


        #Set the velocity and acceleration of the bus.
        self._current_velocity = vf
        self._current_acceleration = a_bus

        # return the power and elapsed time.
        return power_calc, dt
    
    
    def generate_power_frame(self, dist, ext_a):
        '''
        generate_power_frame() takes in distance and external accelerations due to
        hills and friction, and generates an acceleration profile that is limited by max power.
        
        Parameters:
        dist: distance to travel, in meters
        ext_a: external acceleration due to hills/gravity and friction
        
        Returns:
        pandas dataframe of acceleration profile that is limited to travel distance and 
        maximum power. 
        '''
        # Get the raw acceleration profile and its length
        acc_prof = self.get_accel_profile().copy()
        profile_length = len(acc_prof)
        # get the current mass
        m = self._current_mass
        
        # rename acceleration column and time column
        acc_prof=acc_prof.rename(columns = {'accel.[m/s^2]':'true_accel'})
        acc_prof['t_i'] = acc_prof['time[s]']
        del(acc_prof['time[s]'])
        
        # Create change in time column
        acc_prof['t_f'] = acc_prof['t_i'].shift(-1)
        acc_prof['dt'] = (acc_prof['t_f'] - acc_prof['t_i']) 
        del(acc_prof['t_i'])
        del(acc_prof['t_f'])
        
        # create columns for initial velocity, power, and acceleration
        acc_prof['vi'] = 0.0
        acc_prof['power'] = 0.0
        acc_prof['a_param'] =0.0

        
        # create varables for travel distance and iteratable for the while loop
        trav_dist = 0
        i = 0
                    
        # Calculate the maximum power    
        pmax = self._max_power * 1000
        # while the iteratable is less than the profile length,
        while (i<profile_length):
            # get the current true acceleration, change in time, and velocity
            tr_acc = acc_prof['true_accel'][i]
            dt = acc_prof['dt'][i]
            vi = acc_prof['vi'][i]
            
            # add the wind acceleration to the external accelerations:
            net_outside_a = ext_a + self.get_aerodynamic_drag(vi, 0, 1.2)
            
            # do a preliminary power calculation
            P = self.get_power(tr_acc, dt, vi, m, net_outside_a, self._i_factor)
            
            # Chack if the power exceeds the limit:
            if P > (pmax):
                # recalculate the maximum acceleration based on the power output maxing out
                tr_acc = self.get_max_accel(pmax, dt, vi, m, net_outside_a, self._i_factor)
                
                # recalculate the power based on the max acceleration
                P = self.get_power(tr_acc, dt, vi, m, net_outside_a, self._i_factor)
                
            # set the column varaibles to the adjustments
            vi += tr_acc * dt
            dx = vi*dt
            
            # adjust the travel distance
            trav_dist += dx

            # Set the next initial velocity to the current velocity
            acc_prof.at[i+1,'vi'] = vi
            
            # Set the true acceleration, power, change in time, and acceleration parameter
            acc_prof.at[i, 'true_accel'] = tr_acc
            acc_prof.at[i, 'power']= P
            acc_prof.at[i, 'dt'] = dt
            acc_prof.at[i, 'a_param'] = (net_outside_a)
            
            # iterate
            i+=1
            
        # drop the tail indicies of the dataframe, return the profile
        acc_prof.drop(acc_prof.tail(2).index,inplace=True)
        return acc_prof
    
    
    def get_power(self, a_true, dt, vi, m, a_out, i_factor):
        '''
        get_power calculates the power it takes to accelerate based on time,
        mass, external acceleration, initial velocity, and intertial factor.
        This probably should be a private method.
        '''
        P = m * (a_out + i_factor * a_true) * ((2 * vi + dt * a_true) / (2))
        return P


    def get_max_accel(self, pmax, dt, vi, m, a_out, i_factor):
        '''
        get_max_accel backcalculates from a maxumum power to the highest acceleration for those conditions.
        This probably should be a private method.
        '''
        sub1 = (4*(i_factor*dt)*((2*pmax/m) - 2*a_out*vi) + (-a_out * dt - 2*vi*i_factor)**2)
        num = (sub1)**(1/2) - (a_out * dt) - (2 * vi * i_factor)
        denom = 2 * i_factor * dt

        acc = num/denom

        if acc < 0:
            acc = -acc
    
        return acc
    
    
    def accelerate_v5(self, dist, ext_a):
        '''
        accelerate_v5() takes in distance traveled and the external 
        accelerational forces, and calculates the power needed to
        accelerate the bus for that distance, as well as the elapsed time.
        
        Parameters:
        dist: distance of travel in meters as a float.
        ext_acc: external accelerative force the bus is currently experiencing.
        
        Returns:
        a touple containing the power [W] used and the elapsed time [s]. 
        '''
        
        # Generate the acceleration profile for this set of conditions and the current state of the bus
        power_frame = self.generate_power_frame(dist, ext_a)
        
        # Get the length of the acceleration profile
        power_frame_len = len(power_frame)
        
        # Get the last velocity in the power frame as the maximum velocity
        profile_max = list(power_frame['vi'])[-1]
        
        # Get the bus current mass and velocity
        m = self._current_mass
        cur_v = self._current_velocity
        
        # Set up dummy variables to hold the travel distance, powers, and time changes
        trav_dist = 0
        p_list = []
        dt_list = []
        
        # set up an iteratable for the while loop.
        i=0
        
        # set up dimmy variables for p, dt, a, outside_a:
        p, dt, a, outside_a = None, None, None, None
        
        #while the distance to be traveled is less than the traveled distance
        while (dist > trav_dist):
            
            #Check if the velocity is less than the max velocity
            if (cur_v < profile_max):
                
                # if the current velocity is within the acceleration profile,
                # get the index of the closest velocity in the profile
                v_closest_index = list((power_frame['vi']-cur_v).abs().argsort())[0]
                
                # get the power, time, acceleration, and external acceleration for that entry
                p = power_frame['power'][v_closest_index]
                dt = power_frame['dt'][v_closest_index]
                a = power_frame['true_accel'][v_closest_index]
                outside_a = power_frame['a_param'][v_closest_index]
                
            # if the velocity profile is greater than the max velocity:
            else:
                
                # Get the last index of power, dt, and acceleration
                p = list(power_frame['power'])[-2]
                dt = list(power_frame['dt'])[-2]
                a = list(power_frame['true_accel'])[-1]
                
                # add the wind acceleration to the external acceleration
                outside_a = ext_a + self.get_aerodynamic_drag(cur_v, 0, 1.2)

            # Get the dx, dv
            dv = dt*a
            dx = cur_v*dt

            # Check if this would go over the distance
            if dx > abs(trav_dist - dist):

                # get the change in distance to the end of the segment
                dx = abs(trav_dist - dist)

                # get the new change in time
                dt = dx * (self._i_factor * a + outside_a) * m / p

                # calculate the new dv
                dv = a * dt 

            # adjust the current velocity and travel distance
            cur_v += dv
            trav_dist += dx

            #append the power and change in time
            p_list.append(p)
            dt_list.append(dt)
                
            # iterate by 1
            i+=1
        
        # covert to power, time, to series and sum the times
        pow_ser = pd.Series(p_list)
        tee_ser = pd.Series(dt_list)
        travel_time = tee_ser.sum()
        
        # calculate the mean power for the step
        mean_power = (pow_ser * tee_ser).sum() / travel_time
        
        self._current_velocity = cur_v
        
        #return the mean power, travel time
        return mean_power, travel_time
    
    
    def maintain_v5(self, dist, ext_acc):
        '''
        maintain_v5() takes in distance traveled and the external 
        accelerational forces, and calculates the power needed to
        maintain the current velocity and the elapsed time.
        
        Parameters:
        dist: distance of travel in meters as a float.
        ext_acc: external accelerative force the bus is currently experiencing.
        
        Returns:
        a touple containing the power [W] used and the elapsed time [s]. 
        '''
        # get the current velocity and mass
        v_i = self._current_velocity
        m = self._current_mass
        
        # Get the aero drag
        wind_a = self.get_aerodynamic_drag(v_i, 0, 1.2)
        
        # get the maximum power in watts
        max_p =  (self._max_power * 1000)
        
        # set up a dummy variable for net acceleration of the bus
        bus_a = 0
        
        # combine the external accelerations
        a_subsum = ext_acc + wind_a
        
        # Get the preliminary time change
        dt = dist/v_i
        
        # get the initial pass of final velocity
        v_f = v_i
            
   
        # first pass calculation of power
        power_calc = (self._i_factor * bus_a + a_subsum) * m * v_i
        
        # check that the bus brakes are capable of slowing adequately
        if (a_subsum + bus_a) < self._a_braking:
            while (a_subsum+ bus_a) < self._a_braking:
                
                bus_a += .01
                
                # recalculate time based on kinematic equations
                dt = ((2*bus_a*dist + v_i**2)**(1/2) - v_i) / bus_a 
                
                # get the new power calculation
                power_calc = self.get_power(bus_a, dt, v_i, m, a_subsum, self._i_factor) 
                
            
    
        # if the preliminary power calculation exceeds the maximum power
        if power_calc > max_p:
            
            # cycle until the power is matching the maximum power
            while power_calc > max_p:
                
                # decreace the bus' acceleration
                bus_a -= 10/max_p
                
                # recalculate time based on kinematic equations
                dt = ((2*bus_a*dist + v_i**2)**(1/2) - v_i) / bus_a 
                
                # get the new power calculation
                power_calc = self.get_power(bus_a, dt, v_i, m, a_subsum, self._i_factor)     
    
                
                
        
        
        # adjust the current velocity according to the bus acceleration
        self._current_velocity = v_i + bus_a * dt
        
        return power_calc, dt
    
    
    # Getter methods, probably redundant and should turn into just public variables
    
    def get_accel_profile(self):
        return self._accel_profile_df
    
    def velocity(self):
        return self._current_velocity
    
    def set_v(self, v):
        self._current_velocity = v
    
    def max_velocity(self):
        return self._max_v
    
    def update_mass(self, val):
        self._current_mass = val
        return None
    
    def get_n_riders(self):
        return self._passengers
    
    def get_mass(self):
        return self._current_mass
    
    def get_acceleration(self):
        return self._current_accel
    
    def get_fric_coeff(self):
        return self._fric_coeff

    def get_b_accel(self):
        return self._a_braking
    
    def get_motor_eff(self):
        return self._motor_eff
    
    def get_invert_eff(self):
        return self._invert_eff
    
    def get_auxill_efficiency(self):
        return self._eff_aux
    
    def get_aux_load(self):
        return self._aux_load
    
    def get_regen_eff(self):
        return self._regen
        
    