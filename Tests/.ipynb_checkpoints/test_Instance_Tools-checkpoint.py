'''
test_Instance_Tools.py

tests the methods in Instance_Tools.py
'''

import unittest
import sys
sys.path.append('../src/reRoute_Dynamics_Core')
import Instance_Tools as it
import numpy as np

class TestGenerateRiders(unittest.TestCase):
    # Has no handling for negative riderhsip. 
    def test_1_stop(self):
        result = it.generate_riders(1, 5, 1)
        expected = [0]
        self.assertEqual(expected, result)
        
    def test_10_stop(self):
        vals = []
        for i in range(1000):
            result = np.asarray(it.generate_riders(10, 5))
            result = result[result > 0]
            if len(result)>0:
                result = result.mean()
            else:
                result = 0
            vals.append(result)
        mean = np.asarray(vals).mean()
        expected = 1.3

        self.assertEqual(expected, round(mean,1))
        
class TestCheckHitSignal(unittest.TestCase):
    
    def test_hit_signal(self):
        expected = .541666
        results = 0
        for i in range(10000):
            result = it.check_hit_signal(expected, i)
            results += result
        mean = results/10000
        self.assertEqual(0.5458, mean)
            

class TestDetermineStopType(unittest.TestCase):
    
    def test_stop_type(self):
        rider_changes = [0,1,0,-1,0,0,0,1,0,0,0,-1]
        hit_signals =   [1,0,0,0,0,0,0,0,0,0,0,1]
        hit_signs =     [0,0,0,0,0,0,1,0,0,0,0,0]
        
        results = it.determine_stop_type(rider_changes, hit_signals, hit_signs)
        
        self.assertEqual(len(hit_signs), len(results))
        self.assertEqual([1, 1, 0, 1], results[-1])
        
class TestGetStopType(unittest.TestCase):
    
    def test_get_stop_type(self):
        self.assertEqual([0,1,1,1],it.get_stop_type(0,4,'test',-2))
        

class TestGetDistancesToStops(unittest.TestCase):
    
    def test_works(self):
        rider_changes = [0,1,0,-1,0,0,0,1,0,0,0,-1]
        hit_signals =   [1,0,0,0,0,0,0,0,0,0,0,1]
        hit_signs =     [0,0,0,0,0,0,1,0,0,0,0,0]
        dists =         [1,2,3,4,5,6,7,8,9,10,11,12]
        expected =      [0,0,1, 0, 2, 1, 0, 0, 3, 2, 1, 0]
        
        types = it.determine_stop_type(rider_changes, hit_signals, hit_signs)
        results = it.get_distances_to_stops(types, dists)
        self.assertEqual(expected, results)

if __name__ == '__main__':
    unittest.main()