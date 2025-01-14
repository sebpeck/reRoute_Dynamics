'''
test_Physics_Engine.py

contains the tests for Physics_Engine.py

Notes:
1/10/2025 - Maintain needs to get checked for why it could possibly return a v_f of nan. 
'''
import unittest
import sys
import os
sys.path.append('../src/reRoute_Dynamics_Core')
import pandas as pd
import numpy as np
import itertools

import Physics_Engine as pe

# set up the min, max of each parameter
bearings = (0, 360) # deg 0
velocitys = (0, 40) #m/s 1
wind_bearings = bearings # deg 2
wind_speeds = (0, 22) # m/s 3
air_densities = (.9, 1.5) # kg/m^3 4
drag_coefs = (.2, 2) # 5
a_frontals = (.09, 1.5) #m^2, smart car and semi truck 6
grades = (-17, 17) # %, based on seattle 7
masses = (700, 37000) # kg, smart car and loaded semi 8
fric_coefs = (.02, .5) # engineering toolbox 9
a_brakings = (.098, 9.81) # m/s^2, grav 10
br_factors = (0, 1)  # 11
inert_factors = (0.001, 0.0730) # wikipedia 12
max_dist = (300, 1000) #m, estimates 13
dx_travel = (1, 300) #m, estimates 14
max_power = (3000, 60000) # watts, est based on ev smart car, ev semi 15
max_acc = (0.098, 9.81) # m/s^2, grav 16
max_timestep = (.01, 60) # seconds 17
a_prof_path = '../Examples/KC_Example_Data/Acceleration_Profiles/Braunschweig_Acceleration.csv'
a_prof = pd.read_csv(a_prof_path, header=None)
a_prof[1] = a_prof[1]*9.81
a_prof = (a_prof,)


range_list = [bearings, velocitys, wind_bearings, wind_speeds, air_densities, drag_coefs, a_frontals, grades, masses, fric_coefs, a_brakings, br_factors, inert_factors, max_dist, dx_travel, max_power, max_acc, max_timestep, a_prof]

print("Prepping samples...")
param_samples = []
for param_options in list(itertools.product(*range_list)):
    wind_params = param_options[0:7]
    grade_params = param_options[7:10]
    wind_force = pe.calculate_wind_force(*wind_params)
    grade_force = pe.calculate_grade_force(*grade_params)
    bdist_params = (param_options[1], param_options[8], grade_force, wind_force) + param_options[10:14]
    brake_params = (param_options[1], param_options[8], param_options[14], grade_force, wind_force) + param_options[10:14]
    main_params = (param_options[1], param_options[8], param_options[14], grade_force, wind_force) + param_options[10:13] + (param_options[15],)
    acc_params = (param_options[1], param_options[8], param_options[14], grade_force, wind_force,param_options[18]) + param_options[10:13] + param_options[15:18]
    param_samples.append([wind_params, grade_params, bdist_params, brake_params, main_params, acc_params])
print('Samples prepped.')
    
class TestCalculateWindForce(unittest.TestCase):
    print('windforcetests')
    # Legitemately I don't know how to properly test-case this one. 
    # Reasonably, there should be errors resulting for things like:
    # - Negative Frontal area
    # - Negative air density
    # - Negative drag coefficient
    # but so long as the bearing and speed stuff makes sense and the math maths,
    # it should be fine for now. 
    
    def test_no_wind(self):
        # Wind matching angle, no wind speed
        called = pe.calculate_wind_force(0, 10, 0, 0, 1.2, 1, 5)
        expected = 300
        self.assertEqual(called, expected)
        
        
    def test_wind_0(self):
        # Wind matching angle, matching wind speed
        called = pe.calculate_wind_force(0, 10, 0, 10, 1.2, 1, 5)
        expected = -1200
        self.assertEqual(called, expected)
        
        
    def test_wind_180(self):
        # Wind opposite angle, matching speed
        called = pe.calculate_wind_force(0, 10, 180, 10, 1.2, 1, 5)
        expected = 0
        self.assertEqual(called, expected)
        
        
    def test_wind_90(self):
        # Wind side angle, matching speed
        called = pe.calculate_wind_force(0, 10, 90, 10, 1.2, 1, 5)
        expected = -300
        self.assertEqual(called, expected)
        
        
    def test_wind_45(self):
        # Wind 45 angle, matching speed, bus should accel
        called = round(pe.calculate_wind_force(0, 10, 45, 10, 1.2, 1, 5),2)
        expected = -874.26
        self.assertEqual(called, expected)
        
    def test_negative_wind_bearing(self):
        # should give same answer as Wind 45 angle, matching speed, bus should accel
        called = round(pe.calculate_wind_force(0, 10, -315, 10, 1.2, 1, 5),2)
        expected = -874.26
        self.assertEqual(called, expected)
        
    def test_full_loop_negative_wind_bearing(self):
        # should give same answer as wind 45 angle
        called = round(pe.calculate_wind_force(0, 10, 45-360-360, 10, 1.2, 1, 5),2)
        expected = -874.26
        self.assertEqual(called, expected)
        
    def test_bus_moving_backwards(self):
        # This should not ever happen, but here's a test that if it does, the math is still correct.
        # and that the force balances out. 
        called = pe.calculate_wind_force(0, -10, 0, 10, 1.2, 1, 5)
        expected = -0
        self.assertEqual(called, expected)
    
    
class TestSign(unittest.TestCase):
    print('signtests')
    
    def test_negative(self):
        called = pe.sign(-500)
        expected = -1
        self.assertEqual(called, expected)
        
        
    def test_positive(self):
        called = pe.sign(500)
        expected = 1
        self.assertEqual(called, expected)
        
        
    def test_zero(self):
        called = pe.sign(0)
        expected = 1
        self.assertEqual(called, expected)
        

class TestCalculateGradeForce(unittest.TestCase):
    print('gradeforcetests')
    
    # Like above, this should have checks in place to make sure there arent:
    # - Negative or zero mass
    # - Negative or zero friction coefficient
    
    mass = 100
    
    def test_zero_grade(self):
        called = pe.calculate_grade_force(0, self.mass)
        expected = 9.81
        self.assertEqual(called, expected)
        
    
    def test_up_grade(self):
        called = round(pe.calculate_grade_force(10, self.mass),2)
        expected = 107.37
        self.assertEqual(called, expected)
        
        
    def test_down_grade(self):
        called = round(pe.calculate_grade_force(-10, self.mass),2)
        expected = -87.85
        self.assertEqual(called, expected)
    
    def test_no_mass(self):
        called = pe.calculate_grade_force(0, 0)
        expected = 0
        self.assertEqual(called, expected)
        
    def test_negative_mass(self):
        # Not physically possible, but y'never know...
        called = pe.calculate_grade_force(0, -100)
        expected = -9.81
        self.assertEqual(called, expected)
    

class TestGetBrakingDistance(unittest.TestCase):
    print('brakingdistests')
    
    # This should probably have checks in the method
    # to verify that:
    # - velocity isn't negative
    # - mass isn't 0 or below
    # - max distance isn't negative or zero
    # - inertial factor isn't below 1
    
    
    
    mass = 100
    fi = 1.1
    a_br = .40
    fg = 9.81 # Assuming flat, no wind
    fw = 0
    dx = 300
    def test_basic_calc(self):
        # This should not trigger any of the loops and result in the
        # same outcome as a bog-standard kinematic equation.
        result = pe.get_braking_distance(10,
                                         self.mass,
                                         self.fg,
                                         self.fw,
                                         self.a_br,
                                         .5,
                                         self.fi,
                                         self.dx)
        dist = result['dx']
        factor = result['bf']
        ad = result['ad']
        
        expected = self.a_br * .5 * self.fi + (self.fg + self.fw)/self.mass
        expected = ((-10)**2)/(2*expected)
        self.assertEqual(dist, expected)
    
    
    def test_overbraking(self):
        # this SHOULD trigger the first loop, resulting in a decrease in
        # braking factor, resulting in a modulated acceleration
        result = pe.get_braking_distance(10,
                                 self.mass,
                                 self.fg,
                                 self.fw,
                                 self.a_br,
                                 1,
                                 self.fi,
                                 self.dx)
        dist = result['dx']
        factor = result['bf']
        ad = result['ad']
        
        called = round(ad,2)
        expected = round(self.a_br*self.fi, 2)
        self.assertEqual(called, expected)
        
    def test_underbraking_and_neg_bf(self):
        # This should raise the braking factor until the
        # minimum distance dx (in this case 300) is hit.
        result = pe.get_braking_distance(10,
                                 self.mass,
                                 self.fg,
                                 self.fw,
                                 self.a_br,
                                 -1,
                                 self.fi,
                                 self.dx)
        dist = result['dx']
        factor = result['bf']
        ad = result['ad']
        
        called = round(dist,0)
        expected = round(self.dx, 0)
        self.assertEqual(called, expected)
        
    def test_correct_decell_calc(self):
        # this should yeild the same decelleration value,
        # and also check that ad is correctly returned
        result = pe.get_braking_distance(10,
                                 self.mass,
                                 self.fg,
                                 self.fw,
                                 self.a_br,
                                 .5,
                                 self.fi,
                                 self.dx)
        dist = result['dx']
        factor = result['bf']
        ad = result['ad']
        
        expected = self.a_br * .5 * self.fi + (self.fg + self.fw)/self.mass
        self.assertEqual(ad, expected)
        
        
# ---------------------------------------------- General test methods for braking, maintaining, and accelerating
def test_vf_valid(self, vf, params):
    against_min = int(vf<velocitys[0])
    is_valid = int(np.isnan(vf))
    self.assertEqual(is_valid, 0, "vf: {}, ".format(vf)+str(params))
    self.assertEqual(against_min, 0, "vf: {}, ".format(vf)+str(params))
    
def test_t_valid(self, t, params):
    is_possible = int(t<0)
    is_valid = int(np.isnan(t))
    self.assertEqual(is_possible, 0, "t: {}, ".format(t)+str(params))
    self.assertEqual(is_valid, 0, "t: {}, ".format(t)+str(params))
    
def test_P_valid(self, P, params):
    is_valid = int(np.isnan(P))
    self.assertEqual(is_valid, 0, "P: {}, ".format(P)+str(params))
    

class TestBrake(unittest.TestCase):
    print('braketest')
    vi = 10
    mass = 100
    dx = 100
    fg = 9.81 # assuming flat, no wind
    fw = 0
    a_br = .4
    fb = .5
    fi = 1.1
    dx_max = 300
    
    def test_full_stop(self):
        # this should come to a full stop, uses braking distance check
        # to calculate braking distance used
        dist = pe.get_braking_distance(self.vi,
                                self.mass,
                                self.fg,
                                self.fw,
                                self.a_br,
                                self.fb,
                                self.fi,
                                self.dx_max)
        result = pe.brake(self.vi,
                          self.mass,
                          dist['dx'],
                          self.fg,
                          self.fw,
                          self.a_br,
                          dist['bf'],
                          self.fi,
                          self.dx_max)
        
        expected = 0
        
        self.assertEqual(result['v_f'], expected)
        
    def test_no_stop(self):
        # tapping the brake for a distance of zero should mean no velocity change.
        dist = pe.get_braking_distance(self.vi,
                                self.mass,
                                self.fg,
                                self.fw,
                                self.a_br,
                                self.fb,
                                self.fi,
                                self.dx_max)
        result = pe.brake(self.vi,
                          self.mass,
                          0,
                          self.fg,
                          self.fw,
                          self.a_br,
                          dist['bf'],
                          self.fi,
                          self.dx_max)
        
        expected = 0
        self.assertEqual(result['v_f'], self.vi)
    
    
    def test_all_valid(self):
        # This iterates through the combos of min and max possible values for each param and
        # checks if the output are valid responses. 
        for params in param_samples:
            result = pe.brake(*params[3])
            test_vf_valid(self, result['v_f'], params[3])
            test_t_valid(self, result['dt'], params[3])
            test_P_valid(self, result['P'], params[3])

            
class TestMaintain(unittest.TestCase):
    print('maintest')
    def test_all_valid(self):
        # Check if any of the responses are invalid, a la negative time, negative velocity, or nan
        for params in param_samples:
            result = pe.maintain(*params[4])
            test_vf_valid(self, result['v_f'], params[4])
            test_t_valid(self, result['dt'], params[4])
            test_P_valid(self, result['P'], params[4])

    

class TestAccelerate(unittest.TestCase):
    print('acceltest')
    def test_all_valid(self):
        # Check if any of the responses are invalid, a la negative time, negative velocity, or nan
        for params in param_samples:
            result = pe.accelerate(*params[5])
            test_vf_valid(self, result['v_f'], params[5])
            test_t_valid(self, result['dt'], params[5])
            test_P_valid(self, result['P'], params[5])
            
if __name__ == '__main__':
    unittest.main()