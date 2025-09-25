'''
Object_Params.py
S. Peck

Object_Params is used to create, save, and load object classes that store relevant parameters for modeling.

Methods:
load_bus_params() - method to load a Bus object from a saved txt file.
load_ESS_params() - method to load an ESS object from a saved txt file.
load_trip_params() - method ot load a trip object from a saved txt file.
a_eqn() - method to calculate the acceleration at a given time in accordance to a fitting equation.
generate_a_profile() - method to create and save an acceleration profile based on a fitting equation.

Classes:
Bus - a class that is used to store modeling parameters for a bus vehicle. 
copy() - method to create a copy of the bus class
save() - method to save the bus object to a txt file
ESS - a class that is used to store modeling parameters and methods for an Energy Storage System.
copy() - method to create a copy of the ESS
save() - method to save the ESS object to a txt file
bus_E_cap() - method to calculate the energy capacity of the ESS
R_bus() - method to calculate the resistance of the ESS
calc_instance_power() - method to calculate the load on the ESS based on the load needed
calc_voltage_simple() - method to calculate the pack voltage using a simple resistance model at a given power.
Trip - a class that is used to store modeling parameters for a given vehicle trip.
copy() - method to create a copy of the Trip
save() - method to save the trip object to a txt file. 
'''
import pandas as pd
import numpy as np
from ast import literal_eval

import os
import sys
sys.path.insert(0, os.path.abspath("../src/reRoute_Dynamics/"))
A_PROF_PATH = os.path.abspath('../Examples/KC_Example_Data/Acceleration_Profiles/Braunschweig_Acceleration.csv')

class Bus:
    def __init__(self,
                 # Bus Params
                 bus_mass = 13300, # kg, @NF_Excelsior
                 frontal_width = 2.59, # m, @NF_Excelsior
                 frontal_height = 3.38, #m, @NF_Excelsior
                 drag_coeff =.6, # From Aggregate, @Abdelaty_Mohamed
                 friction_coeff = .01, # @Abdelaty_Mohamed
                 braking_accel = 1.5, #m/s^2 @APTA_Braking_Standards <-- Possible source, handbrake road minimum is ~1.5
                 br_factor = .5, # Not tied to any true value.
                 a_factor=.5, # Not tied to any true value
                 i_factor = 1.1, # intertial factor, @Asamer_Et_Al
                 max_dist = 304.8, # m, Assembled from google measurements of offramps from I5
                 a_prof_path = A_PROF_PATH, #@NREL_Drive_Cycles
                 max_acc = .4, # m/s^2, the default accel after profile finishes.
                 max_dt = 1, # s, timestep for the default acceleration betond the profile.
                 max_P = 160000 # Watts, max power output by the motors. Default not tied to a true value.
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
        self.a_max = max_acc
        self.dt_max = max_dt
        self.P_max = max_P
        self.f_a = a_factor
        
            
        self.a_prof[1] = self.a_prof.apply(lambda x: x[1]*9.81, axis=1)
        
        if (self.a_prof.iloc[-1][0] != self.dt_max + self.a_prof.iloc[-2][0]) and self.a_prof.iloc[-1][1] != self.a_max:
            self.a_prof = pd.concat([self.a_prof,  pd.DataFrame([{0: self.dt_max + self.a_prof.iloc[-1][0], 1: self.a_max}])]).reset_index(drop=True)

        
    def copy(self):
        return Bus(self.mass, self._w, self._h, self.Cd, self.Cf, self.a_br, self.f_br, self.f_a, self.f_i, self.dmax, self._a_prof_path, self.a_max, self.dt_max, self.P_max)
    
    def save(self, filepath):
        data = "{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(self.mass,
                                                            self._w,
                                                            self._h,
                                                            self.Cd,
                                                            self.Cf,
                                                            self.a_br,
                                                            self.f_br,
                                                            self.f_a,
                                                            self.f_i,
                                                            self.dmax,
                                                            self._a_prof_path,
                                                            self.a_max,
                                                            self.dt_max,
                                                            self.P_max)
        


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
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8,9, 11, 12, 13]

    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
    args = tuple(data_list)
    
    return Bus(*args)
    

class ESS: 
    def __init__(self,
                 # ESS instance
                 motor_eff = .916, #Approximation
                 inverter_eff = .971, #Approximation
                 aux_eff = .89, # Approximation
                 simple_load = 7000, # Watts - @Abdelaty_Mohamed
                 regen_eff = .6, # @Gallet
                 max_regen = -100000, # Watts,Approximation
                 cell_ocv = 3.3, # Nominal LFP voltage for A123 26650 @A123_26650
                 cell_res = .008, #ohms, cell internal resistance LFP A123 26650 @A123_25550
                 module_struct = (12, 8),# Series, parallel config for a module of cells
                 bus_struct = (16, 1), # Series, parallel config for the modules of the bus
                 cell_cap = 2.3 # Ah, nominal cell capacity @A123_26650
                ):
        self.Em = motor_eff
        self.Ei = inverter_eff
        self.Ea = aux_eff
        self.Er = regen_eff
        self.P_aux = simple_load
        self.P_regen = max_regen
        self.V_cell = cell_ocv
        self.R_cell = cell_res
        self.module_S_P = module_struct
        self.bus_S_P = bus_struct
        self.Q_cell = cell_cap
    
    def bus_E_cap(self):
        return self.V_cell*self.module_S_P[0]*self.bus_S_P[0]*self.Q_cell*self.module_S_P[1]*self.bus_S_P[1]
        
        
    def R_bus(self):
        return self.R_cell*self.module_S_P[0]/self.module_S_P[1]*self.bus_S_P[0]/self.bus_S_P[1] 
    
    def calc_instance_power(self, value):
        '''calc_instance_power takes in a power value,
        and converts it to the corresponding load on the ESS.
        This is a simple stopgap.

        :param value: a power value in Watts, as an int or float.

        :return: converted battery power as a float.
        '''

        # set the battery power to zero.
        bat_pow = 0
        # Including Auxilliary load, though not strictly important at the moment. 
        if (value >= 0):
            # Discharging, converting the needed power into power battery must exert
            bat_pow = value/(self.Em*self.Ei) + (self.P_aux/self.Ea)
        elif(value*self.Er*self.Em > self.P_regen):

            #charging, the regenerative braking ALL the time, max regen is 100
            bat_pow = (value*self.Er*self.Em) + (self.P_aux/self.Ea)

        else:
            bat_pow = self.P_regen + (self.P_aux/self.Ea)


        # Return the battery power.
        return bat_pow
    
    
    def calc_voltage_simple(self, value):
        '''
        Use a simple resistance model to calculate the voltage of a cell based off of a given power.
        '''

        if value <0:
            sign = -1
        else:
            sign = 1
            

        I = sign*np.sqrt(abs(value)/self.R_bus())
        v_module =(self.V_cell - (self.R_cell * I)/(self.bus_S_P[1] * self.module_S_P[1]))
        return v_module
    
    def _calc_voltage_BV(self, value):
        '''
        Use the Butler-Volmer equation to predict cell voltage from a given power, assuming T=25 C and symmetric cell
        Notes:
        9/3/2025 - beta constant needs to be a variable, including changing temperature. This should always assume a symmetric li-ion cell for the math to work.
                   Also needs ot be tested. 
        '''

        # constant for symmetric cell at 25 C
        beta_const = .5*1* 9.648*10**4 /( 8.314 * 25 + 273.15 )
        p_cell = value/self.bus_s_p[1]*self.module_s_p[1]*self.bus_s_p[0]*self.module_s_p[0]
        
        v_cell = np.arcsinh(np.sqrt(p_cell/self.R_cell)/(self.V_cell/self.R_cell*2))/beta_const - self.V_cell
        
        return v_cell
        
        
    def decay_by_c_rate(self, c_rate):
        '''
        Determine the capacity decay coefficient from the c-rate.
        '''
        c=abs(c_rate)
        '''


        if c<5:
            return .01*c/100
        elif c>5:
            return .15*c/100
        else:
            return 0
            
        if c < .1:
            return 0
        elif .1 <= c < .5:
            return .000011*c
        elif .5 <= c < 1.5:
            return .000022*c
        else:
            return .000066*c
        '''
        
        # Per Dan's Advisment, I should just use a linear function for this, rather than a piecewise
        #.000042
        return (c)*.00005


            
        
        
    def copy(self):
        return ESS(self.Em,
                   self.Ei,
                   self.Ea,
                   self.P_aux,
                   self.Er,
                   self.P_regen,
                   self.V_cell,
                   self.R_cell,
                   self.module_S_P,
                   self.bus_S_P,
                   self.Q_cell)
    
    
    def save(self, filepath):
        data = "{}|{}|{}|{}|{}|{}|{}|{}|{}|{}|{}".format(self.Em,
                                             self.Ei,
                                             self.Ea,
                                             self.P_aux,
                                             self.Er,
                                             self.P_regen,
                                             self.V_cell,
                                             self.R_cell,
                                             self.module_S_P,
                                             self.bus_S_P,
                                             self.Q_cell)
        
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
    data_list = data.split('|')
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 10]
    tuple_indexes = [8, 9]
    

    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
        
    for index in tuple_indexes:
        data_list[index] = literal_eval(data_list[index])
        
    args = tuple(data_list)

    return ESS(*args)
          
        
class Trip:
    def __init__(self,
                 pass_mass = 70, # kg Approximation
                 limit_MOE = 4.47, #m/s, or 10 mph <- Custom
                 signal_rest = 65/2, # s, based on @WSDOT_Signals
                 signal_chance=.541666, # Approximation
                 stop_rest = 7, # seconds per passenger boarding, Approximation
                 sign_rest = 7, # Approximation
                 end_rest = 10, # seconds before engine shutdown/startup - Approxiamaiton
                 air_density = 1.2 , # kg/m^3, approximate at 20 deg C.
                 wind_speed = 1.78, # m/s, typical for seattle @Weatherspark_Seattle
                 wind_heading = 'SE', #per @Weatherspark_Seattle
                 temperature = 12.788, #deg C, @Weatherspark_Seattle
                 interp_length = 10, # meters, Approximation
                 mean_ridership = 3.5, # People. Approximate average value across all trips. Not Accurate. 
                 seed = 42, # random seed
                 lg=43, #savgol param
                 deg = 3, #savgol param
                 stop_margin = 1, #m/s, margin for the velocity to be considered 'stopped'
                 traffic = 0 # fraction from zero to 1 representing how slowed the
                             # bus is due to traffic
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
        self.traffic = traffic
        
    def copy(self):
        return Trip(self.m_pass, self.MOE, self.chance_sig, self.t_sig, self.t_stop, self.t_sign, self.t_end, self.p_air, self.T_air, self.wind_bearing, self.v_wind, self.d_interp, self.m_riders, self.seed, self.lg, self.deg, self.stop_margin, self.traffic)
    
    def save(self, filepath):
        data = "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(self.m_pass,
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
                                                                             self.stop_margin,
                                                                             self.traffic)
        
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
    
    numerical_indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17]
    for index in numerical_indexes:
        data_list[index] = float(data_list[index])
    args = tuple(data_list)
    
    return Trip(*args)
    
    

def a_eqn(t, m=-4.9661, b=2.9465):
    '''a_eqn is used to calculate the acceleration at a given time during the acceleration process from zero.
    the default values are based on a fit of the Braunschweig drive cycle.
    
    :param t: time, in seconds, since the bus began accelerating, as a float
    :param m: slope value of the linear fit of 1/t vs ln(v) using data aggregated from Braunschweig https://www.nrel.gov/transportation/drive-cycle-tool/
        Default of -4.9661.
    :param b: intercept value of aformentioned fit as float.
        Default of 2.9465.
    
    :return: acceleration in m/s^2.
    '''
    if t <= 0:
        return 0
    else:
        return -m/(t**2) * np.exp((m/t) + b)

def generate_a_profile(filepath, m=-4.9661, b=2.9465, start=0, stop=34, step=.5):
    '''generate_a_profile() takes the fit parameters for an acceleration profile, and
    generates one for a given range and step and saves at a filepath.
    
    :param filepath: savefile path and filename.
    :param m: slope value of the linear fit of 1/t vs ln(v) using data aggregated from 
    Braunschweig https://www.nrel.gov/transportation/drive-cycle-tool/
        Default of -4.9661.
    :param b: intercept value of aformentioned fit as float. 
        Default of 2.9465.
    :param start: starting value for range.
        Default value of 0
    :param stop: stop value for range. 
        Default value of 34
    :param step: step size for range. 
        Default of .5.
    
    :return: filepath to generated acceleration profile
    '''
    
    # Get the range of times as a series
    total_times = pd.Series(np.arange(start, stop, step))
    
    # Apply the acceleration fit and then combine into a single dataframe
    a_prof = pd.concat([total_times, total_times.apply(lambda x: a_eqn(x, m, b))], axis=1)
    
    # fix to fit the proper acceleration profile formatting and return
    a_prof[1] = a_prof[1].shift(-1)/9.81
    a_prof = a_prof[:-1]
    a_prof.to_csv(filepath, index=False, header=False)
    return filepath