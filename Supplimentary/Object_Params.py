'''
Object_Params.py
'''
import pandas as pd
import numpy as np


class Bus:
    def __init__(self,
                 # Bus Params
                 bus_mass = 13300, # kg, https://en.wikipedia.org/wiki/New_Flyer_Xcelsior
                 frontal_width = 2.59, # m, https://en.wikipedia.org/wiki/New_Flyer_Xcelsior
                 frontal_height = 3.38, #m, https://en.wikipedia.org/wiki/New_Flyer_Xcelsior
                 drag_coeff =.6, # From Aggregate, Abdelaty & Mohamed
                 friction_coeff = .01, # Abdelaty & Mohamed
                 braking_accel = 1.5, #m/s^2 # https://www.apta.com/wp-content/uploads/APTA-BTS-BC-RP-001-05_Rev1.pdf <-- Possible source, handbrake road minimum is ~1.5
                 br_factor = .5, # Not tied to any true value.
                 i_factor = 1.1, # intertial factor, Asamer Et. Al
                 max_dist = 304.8, # m, Assembled from google measurements of offramps from I5
                 a_prof_path = '/media/sebastian/Slepnir/Route_Data/Accel_Prof/acceleration.csv', #TODO: NEEDS SOURCE
                 max_acc = .4, # m/s^2, the default accel after profile finishes.
                 max_dt = 1, # s, timestep for the default acceleration betond the profile.
                ):
        self._a_prof_path=a_prof_path
        self._w = frontal_width
        self._h = frontal_height
        self.mass = bus_mass
        self.area = frontal_width*frontal_height
        self.Cd = drag_coeff
        self.Cf = friction_coeff
        self.a_br = braking_accel
        self.f_br = br_factor
        self.f_i = i_factor
        self.dmax = max_dist
        self.a_prof = pd.read_csv(a_prof_path, header=None)
        self.a_prof[1] = self.a_prof.apply(lambda x: x[1]*9.81, axis=1)
        self.a_max = max_acc
        self.dt_max = max_dt
        
    def copy(self):
        return Bus(self.mass, self._w, self._h, self.Cd, self.Cf, self.a_br, self.f_br, self.f_i, self.dmax, self._a_prof_path, self.a_max, self.dt_max)
    
    def save(self, filepath):
        data = "{},{},{},{},{},{},{},{},{},{},{},{}".format(self.mass,
                                                            self._w,
                                                            self._h,
                                                            self.Cd,
                                                            self.Cf,
                                                            self.a_br,
                                                            self.f_br,
                                                            self.f_i,
                                                            self.dmax,
                                                            self._a_prof_path,
                                                            self.a_max,
                                                            self.dt_max)
        


        # Clear the file
        open(filepath, 'w').close()
        
        # Write the params
        with open(filepath, 'w') as f:
            f.write(data)
        return filepath
    

def load_bus_params(filepath):
    data = ''
    with open(filepath, 'r') as f:
        data = f.read()
    data_list = data.split(',')
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11]
    
    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
    args = tuple(data_list)
    
    return Bus(*args)
    

class ESS: 
    def __init__(self,
                 # ESS instance
                 motor_eff = .916, # <-- NEEDS SOURCE
                 inverter_eff = .971, # <-- NEEDS SOURCE
                 aux_eff = .89, # <-- NEEDS SOURCE
                 simple_load = 7000, # Watts - Adelaty & Mohamed
                 regen_eff = .6, # Gallet Et Al
                 max_regen = -100000, # Watts, currently a guess. Needs Source.
                 max_power = 160000 #w, peak power200kw. 160kw cont.  HDS200 https://cptdb.ca/wiki/index.php/BAE_Systems#HybriDrive
                ):
        self.Em = motor_eff
        self.Ei = inverter_eff
        self.Ea = aux_eff
        self.Er = regen_eff
        self.P_aux = simple_load
        self.P_regen = max_regen
        self.P_max = max_power
        
    def copy(self):
        return ESS(self.Em, self.Ei, self.Ea, self.Er, self.P_aux, self.P_regen, self.P_max)
    
    
    def save(self, filepath):
        data = "{},{},{},{},{},{},{}".format(self.Em, self.Ei, self.Ea, self.Er, self.P_aux, self.P_regen, self.P_max)
        
        # Clear the file
        open(filepath, 'w').close()
        
        # Write the params
        with open(filepath, 'w') as f:
            f.write(data)
        return filepath
    

def load_ESS_params(filepath):
    data = ''
    with open(filepath, 'r') as f:
        data = f.read()
    data_list = data.split(',')
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6]
    
    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
    args = tuple(data_list)
    
    return ESS(*args)
          
        
class Trip:
    def __init__(self,
                 pass_mass = 70, # kg
                 limit_MOE = 4.47, #m/s, or 10 mph <- Custom
                 signal_rest = 65/2, # s, based on https://wsdot.wa.gov/travel/operations-services/traffic-signals
                 signal_chance=.541666,
                 stop_rest = 7, # seconds per passenger boarding
                 sign_rest = 7, # guess.
                 end_rest = 10, # seconds, doesn't need to be precise.
                 air_density = 1.2 , # kg/m^3, approximate at 20 deg C.
                 wind_speed = 1.78, # m/s, typical for seattle per  https://weatherspark.com/y/913/Average-Weather-in-Seattle-Washington-United-States-Year-Round
                 wind_heading = 'SE', #per https://weatherspark.com/y/913/Average-Weather-in-Seattle-Washington-United-States-Year-Round
                 temperature = 12.788, #deg C, per https://weatherspark.com/y/913/Average-Weather-in-Seattle-Washington-United-States-Year-Round 
                 interp_length = 10, # meters, guess. 
                 mean_ridership = 3.5, # People. Back of napkin. 
                 seed = 42, # random seed
                 lg=43, #savgol param
                 deg = 3, #savgol param
                 stop_margin = 1 #m/s, margin for the velocity to be considered 'stopped'
                ):
        self.m_pass = pass_mass
        self.MOE = limit_MOE
        self.chance_sig = signal_chance
        self.t_sig = signal_rest
        self.t_stop = stop_rest
        self.t_sign = sign_rest
        self.t_end = end_rest
        self.p_air = air_density
        self.T_air = temperature
        self.wind_bearing = wind_heading
        self.v_wind = wind_speed
        self.d_interp = interp_length
        self.m_riders = mean_ridership
        self.seed = seed
        self.lg = lg
        self.deg = deg
        self.stop_margin = stop_margin
        
    def copy(self):
        return Trip(self.m_pass, self.MOE, self.chance_sig, self.t_sig, self.t_stop, self.t_sign, self.t_end, self.p_air, self.T_air, self.wind_bearing, self.v_wind, self.d_interp, self.m_riders, self.seed, self.lg, self.deg, self.stop_margin)
    
    def save(self, filepath):
        data = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(self.m_pass,
                                                                             self.MOE,
                                                                             self.chance_sig,
                                                                             self.t_sig,
                                                                             self.t_stop,
                                                                             self.t_sign,
                                                                             self.t_end,
                                                                             self.p_air,
                                                                             self.T_air,
                                                                             self.wind_bearing,
                                                                             self.v_wind,
                                                                             self.d_interp,
                                                                             self.m_riders,
                                                                             self.seed,
                                                                             self.lg,
                                                                             self.deg,
                                                                             self.stop_margin)
        
        # Clear the file
        open(filepath, 'w').close()
        
        # Write the params
        with open(filepath, 'w') as f:
            f.write(data)
        return filepath
    

def load_trip_params(filepath):
    data = '' 
    with open(filepath, 'r') as f:
        data = f.read()
    data_list = data.split(',')
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16]
    
    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
    args = tuple(data_list)
    
    return Trip(*args)
    