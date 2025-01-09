'''
test_ESS.py

tests the methods in ESS.py

Notes:
1/18/2025 - May want to split the instance power calc to a handful of others. Don't wanna spend too much time on unit tests as ESS is very subject to change at the moment. 
'''

import unittest
import sys
sys.path.append('../src/reRoute_Dynamics_Core')
import ESS


class TestInstancePower(unittest.TestCase):
    
    #should ESS metod be split to multiple individual calculations, or is that less efficient?
    def test_in_range(self):
        called = round(ESS.calc_instance_power(0), 2)
        expected = 7865.17
        self.assertEqual(called, expected)
        
        
    def test_below_range(self):
        called = round(ESS.calc_instance_power(-10000000), 2)
        expected = round(-100000 + (7000)/.89, 2)
        self.assertEqual(called, expected)
    
    
    def test_above_range(self):
        called = round(ESS.calc_instance_power(10000000), 2)
        expected = round(10000000/(.916*.971) + (7000)/.89, 2)
        self.assertEqual(called, expected)
        
    def test_within_charge_limit(self):
        called = round(ESS.calc_instance_power(-50000), 2)
        expected = -19614.83
        self.assertEqual(called, expected)
        
if __name__ == '__main__':
    unittest.main()