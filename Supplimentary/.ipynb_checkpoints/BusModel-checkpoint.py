import pandas as pd
import numpy as np

class BusModel:
    def __init__(self,
                 acceleration_profile, # two column dataframe containing time(s), accel(m/s^2)
                 raw_mass = 13041, # kilograms, unknown bus model
                 bus_width = 2.6, # meters, unknown bus model
                 bus_height = 3.3, # meters, unknown bus model
                 drag_coeff = 0.6, # default from unknown bus model
                 wheel_rad = 0.5, # meters, unknown bus model
                 factor = 1.1, # intertial factor.
                 fric_coeff = .01, # default from erica's
                 motor_eff = .916, # unknown source
                 invert_eff = .971, # unknown source
                 max_power = 160, # kW, output motor
                 regen = .6, # unknown source
                 eff_aux = .89, # unknown source, Auxulliary system efficiency?
                 aux_load = 0, # W, no default system load
                 a_braking = -1.5, #m/s^2
                 final_a = .4, # m/s^2, defualt acceleration after profile finishes
                 max_velocity = 26.8224, # m/s, = 60 mph
                 maintain_acceleration = False, # boolean for if bus should maintain the last a in 
                                                # profile or use final_a for extrapolating new vals
                 num_starting_passengers = 0,
                 pass_ave_mass = 80, #kg
                
                 ):
        # raw bus characteristics -- Not all of these are used, but may be useful for the future
        self._empty_mass = raw_mass
        self._bus_width = bus_width
        self._bus_height = bus_height
        self._bus_front_area = self._bus_width * self._bus_height
        self._drag_coeff = drag_coeff
        self._wheel_rad = wheel_rad
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
    
    
    def get_accel_profile(self):
        return self._accel_profile_df
    
    
    def velocity(self):
        return self._current_velocity
    
    
    def set_v(self, v):
        self._current_velocity = v
    
    
    def max_velocity(self):
        return self._max_v
    
    def accelerate(self, dist):
        '''
        accelerate() is a method of calculating the effect of 
        acceleration has on the bus through the acceleration profile.
        Note: This is currently unused in the Trip class, as it is
        not as accurate to the true experience as accelerate_v2().
        
        Parameters:
        dist: distance, in meters, the bus accelerates
        
        Returns:
        power_used: a value, in watts, of power used to accelerate as a float.
        '''
        # Get current velocity
        prev_v = self._current_velocity
        
        # get a copy of the acceleration profile
        prof = self.get_accel_profile().copy()
        
        prof['force[N]']=self._current_mass*prof['accel.[m/s^2]']
        prof['force[N]']=prof['force[N]'].cumsum()
        
        prof['energy[J]'] = prof['force[N]']*(prof['dist.[m]'].shift(-1))
        prof['power[W]'] = prof['energy[J]']*(prof['dt'].shift(-1))
        
        # get the closest velocity in the acceleration profile
        closest_v = prof.iloc[(prof['vel.[m/s]'] - prev_v).abs().argsort()[:1]]
        
        
        # get the index of the closest velocity in the profile
        starting_index = (list(closest_v.index))[0]
        
        # generate an acceleration distance column based on the starting point of closest v
        prof['accel_dist'] = prof['dist.[m]'] - int(closest_v['dist.[m]']) - dist
        
        # zero the velocity to the current velocity
        zeroed_v = prof.iloc[(prof['accel_dist'] - 0).abs().argsort()[:1]]
        
        # get the final index of the zeroed v
        final_index = (list(zeroed_v.index))[0]
        
        # set the current velocity to the final index of the acceleration
        self._current_velocity = prof['vel.[m/s]'][final_index]
        force_applied = prof['force[N]'][final_index] - prof['force[N]'][starting_index] 
        time_taken = prof['time[s]'][final_index] - prof['time[s]'][starting_index] 
        
        # set the distance traveled to the sum of the previous value and the
        # travel distance
        self._distance_traveled = dist + self._distance_traveled
        
        # change bus status
        self._bus_status = 'moving'
        
        # get the raw acceleration energy
        accelerate_energy_raw = force_applied*dist
        
        # Calculate the power used in watts
        power_used = accelerate_energy_raw/time_taken # Watts
        
        # Return the power used
        return power_used
    
    def accelerate_v2(self, dist):
        '''
        accelerate_v2() takes in a distance value, and determines the power
        and time change the bus requires to accelerate that distance.
        
        Parameters:
        dist: distance, in meters (float), the bus accelerates.
        
        Returns:
        a touple containing the power [W] used and the elapsed time [s]. 
        '''
        # get the current velocity of the bus.
        prev_v = self._current_velocity
        
        # get a copy of the acceleration profile.
        prof = self.get_accel_profile().copy()

        # set the previous velocity as the zero index in velocity
        prof.at[0, 'vel.[m/s]'] = prev_v 
        
        # Set up columns for force, energy, power, and distance change
        prof['force[N]']=self._current_mass*prof['accel.[m/s^2]']
        prof['energy[J]'] = 0
        prof['power[W]'] = 0
        prof['dx[m]'] = 0
        
        # set up variables for accelerated distance and index
        accelerated_dist = 0
        i=0
        
        # Loop through the acceleration profile so long as the bus has not yet traveled the
        # target distance and so long as the power does not exceed the maximum output of the
        # motors. 
        while ((dist > accelerated_dist) & (prof['power[W]'].iloc[i] < self._max_power*1000)):
            # Set v0 to the velocity at index
            v0 = prof['vel.[m/s]'][i]
            
            # get current acceleration
            a = prof['accel.[m/s^2]'][i]
            
            # get the change in time
            dt = prof['dt'][i+1]
            
            # using kinematics equations, calculate change in distance,
            # and the final velocity
            dx = v0*dt + (1/2)*a*(dt**2)
            vf = np.sqrt(v0**2 + 2*a*dx)
            
            # Set the next variables to the newly calculated ones.
            prof.at[i+1, 'vel.[m/s]'] = vf
            prof.at[i+1, 'dx[m]'] = dx
            prof.at[i+1, 'energy[J]'] = prof['force[N]'][i+1]*dx
            prof.at[i+1, 'power[W]']= prof['energy[J]'].iloc[i+1]/dt
            accelerated_dist += dx
            
            # Index by one.
            i+=1
            
        # Get the final velocity value and elapsed time
        final_v = prof['vel.[m/s]'].iloc[i]
        dt = prof['time[s]'].iloc[i]
        
        # set the current values to the bus's storage variables
        self._current_velocity = final_v
        self.distance_traveled = dist + self._distance_traveled
        self._bus_status = 'moving'
        
        # Calculate the mean power used in watts
        power_used = prof['power[W]'].cumsum()[i]/i # Watts
        
        # return the power used and the elapsed time
        return power_used, dt
        
        
    def brake(self, dist, braking_factor):
        '''
        brake() takes in a distance and braking factor,
        and determines the power and elapsed time for the bus to
        experience those conditions.
        '''
        
        # Get the current velocity and braking acceleration, 
        # using the factor to help determine how intense the braking is.
        v0 = self._current_velocity
        braking_acc = self._a_braking*braking_factor
        
        # Use the kinematics equations to calculate the final velocity
        v2 = v0**2
        dv2 = braking_acc*dist*2
        final_v = 0
        if (v2>np.abs(dv2)):
            #print(v2 + dv2)
            final_v = np.sqrt(v0**2 + 2*braking_acc*dist)
        self._current_velocity = final_v
        
        # get the braking force and energy, and set the current acceleration
        # to the braking acceleration. 
        braking_force = braking_acc*self._current_mass
        braking_energy_raw = braking_force*dist
        self._current_accel = braking_acc
        
        # Calculate the power of the brakign and the change in distance
        power_calc = braking_energy_raw/((final_v - v0)/braking_acc)
        dt = dist / (abs(final_v - v0)/2)
        
        # return the power and elapsed time.
        return power_calc, dt
    
    
    def maintain(self, dist, ext_acc):
        '''
        maintain() takes in distance traveled and the external 
        accelerational forces, and calculates the power needed to
        maintain the current velocity and the elapsed time.
        
        Parameters:
        dist: distance of travel in meters as a float.
        ext_acc: external accelerative force the bus is currently experiencing.
        
        Returns:
        a touple containing the power [W] used and the elapsed time [s]. 
        '''
        
        # The current external acceleration from what was passed
        current_ext_acc = ext_acc
        
        ## THIS RIGHT HERE, CHIEF! THIS IS THE PROBLEM!!!!

        # Get current bus mass
        cur_mass = self._current_mass
        
        # get the external force from the mass and acceleration
        ext_force = cur_mass * current_ext_acc
        
        # get the force of inertia
        in_force = self.get_inertial_force()
        
        # The engine force should compensate to balance
        # the difference between the inertia and the external force
        eng_force = in_force+(ext_force)
        
        # get the energy to maintain the force for the distance
        maintain_energy_raw = eng_force*dist
      
        dt = dist/self._current_velocity
        
        power_calc = maintain_energy_raw/dt
        
        # return the energy value
        return power_calc, dt
    
    
    def get_braking_distance(self, velocity, ext_acc):
        '''
        get_braking_distance() takes in a velocity and the current external accelerations,
        and provides a calculation as to how far the bus would have to brake to reach a velocity of zero.
        
        Parameters:
        velocity: velocity, in m/s, the bus is travelling at
        ext_acc: the external accelerative force in m/s^2 the bus is experiencing.
        
        Returns:
        a distance, in meters, of how far the bus would have to brake to reach a velocity of zero.
        '''
        
        # Using the kinematics equation vf^2 = vi^2 + 2a(dX) to get
        # braking distance
        return -1*round((velocity**2 / (2*(self._a_braking - ext_acc))), 5) # Convert to meters

    
    def get_aerodynamic_drag(self, wind_speed, air_density):
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
        air_drag = (self._drag_coeff * self._bus_front_area * (air_density/2) * (self._current_velocity - wind_speed)**2)/self._current_mass
        
        # Return said value.
        return air_drag # acceleration of air drag
    
    
    def get_inertial_force(self):
        '''
        get_inertial_force calculates the current force of intertia the bus is undergoing.
        
        Parameters:
        N/A
        
        Returns:
        inertial force of the bus, in Newtons.
        
        TODO:
        Ask about Intertial Factor - What does it represent? Rolling resistancE?
        '''
        
        # Use the intertial factor, current mass, and current acceleration,
        # to et inertial force.
        inertia = self._i_factor * self._current_mass * self._current_accel
        return inertia # Force of inertia
    
    
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
    
    
    # Getter methods
    def update_mass(self, val):
        self._current_mass = val
        return None
    
    def get_n_riders(self):
        return self._passengers
    
    def get_mass(self):
        return self._current_mass
    
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