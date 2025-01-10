'''
test_Object_Params.py

contains testing classes and cases for Object_Params.py

Notes:
1/8/2025 - This could do with some additional tests, but as it stands OP is simple enough that 
           very little can actually go wrong with it. 
1/9/2025 - I have since discovered something that can go wrong with this. There need to be tests that can validate that parameters do NOT get switched around when saving and loading, as that can be a bit of an issue. 
'''

import unittest
import sys
sys.path.append('../src/reRoute_Dynamics_Core')
import Object_Params as op
import pandas as pd
from pandas.testing import assert_frame_equal


# Testing Bus Class
class TestBus(unittest.TestCase):
    
    
    def test_can_save_load(self):
        bus = op.Bus()
        bus.mass = 999999
        bus.save('test_bus.txt') 
        loaded_bus = op.load_bus_params('test_bus.txt')
        
        self.assertEqual(bus.mass, loaded_bus.mass)
        
        
    def test_can_copy(self):
        bus = op.Bus()
        bus.Cd = 1
        self.assertEqual(bus.copy().Cd, bus.Cd)
        
        
    def test_can_load_a_prof(self):
        bus = op.Bus(a_prof_path = "../Examples/KC_Example_Data/Acceleration_Profiles/Braunschweig_Acceleration.csv")
        self.assertEqual(len(bus.a_prof), 67)
        
        
# Testing the ESS class.
class TestEss(unittest.TestCase):
    
    def test_can_save_load(self):
        ESS = op.ESS()
        ESS.Em = 999999
        ESS.save('test_ESS.txt') 
        loaded_ESS = op.load_ESS_params('test_ESS.txt')
        self.assertEqual(ESS.Em, loaded_ESS.Em)
        
        
    def test_can_copy(self):
        ESS = op.ESS()
        ESS.P_regen= 1
        self.assertEqual(ESS.copy().P_regen, ESS.P_regen)

        
# Testing the Trip class.
class TestTrip(unittest.TestCase):
    
    def test_can_save_load(self):
        trip = op.Trip()
        trip.m_pass = -1
        trip.save('test_trip.txt')
        loaded_trip = op.load_trip_params('test_trip.txt')
        self.assertEqual(trip.m_pass, loaded_trip.m_pass)
        
    def test_can_copy(self):
        trip = op.Trip()
        trip.seed = 4184
        
        self.assertEqual(trip.copy().seed, trip.seed)
    

# Testing the Acceleration Profile Generation.
class TestAProfGen(unittest.TestCase):
    
    
    def test_a_eqn_valid(self):
        t = 1
        m = 1
        b = 1
        res = -7.389
        self.assertEqual(round(op.a_eqn(1, 1, 1),3), -7.389)
        
        
    def test_a_eqn_out_time(self):
        t = -4
        self.assertEqual(op.a_eqn(t), 0)
        
        
    def test_a_profile_braunschweig(self):
        a_prof_path = "../Examples/KC_Example_Data/Acceleration_Profiles/Braunschweig_Acceleration.csv"
        braun = pd.read_csv(a_prof_path, header=None)
        generated = pd.read_csv(op.generate_a_profile('test_profile.csv'), header=None)
        assert_frame_equal(generated, braun)
        
        
        
if __name__ == '__main__':
    unittest.main()
        
        