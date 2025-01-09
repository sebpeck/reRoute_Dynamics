'''
test_Geography_Tools.py

tests the methods in Geography_Tools.py
'''

import unittest
import sys
sys.path.append('../src/reRoute_Dynamics_Core')
import Geography_Tools as gt

class TestHaversineFormula(unittest.TestCase):
    
    def test_large_range(self):
        called = round(gt.haversine_formula(0, 0, 5, 5),2)
        expected = 785.33
        self.assertEqual(called, expected)
        
        
    def test_small_range(self):
        called = round(gt.haversine_formula(47.68688, -122.37682, 47.68678, -122.37682)*1000, 2)
        expected = 11.11
        # Google clocks this as 10.62 meters. Other calculators have it as 11.11 Not sure what to make of that.
        self.assertEqual(called, expected)
    
    
    def test_identical(self):
        # Having this return zero can cause issues, so I want it to make sure it's at least A value. 
        called = gt.haversine_formula(47.68679, -122.37682, 47.68679, -122.37682)
        expected = .0001
        self.assertEqual(called, expected)
        
        
class TestCompassHeading(unittest.TestCase):
    
class TestHeadingToAngle(unittest.TestCase):
    
class TestPointBearing(unittest.TestCase):
    
class TestBoundingBos(unittest.TestCase):
    
class TestInterpolatePoints(unittest.TestCase):
    
class TestRepeatIDRemover(unittest.TestCase):
    
class TestVerboseLineUpdater(unittest.TestCase):
    
class TestQueryElevationChanges(unittest.TestCase):
    
class TestQueryStops(unittest.TestCase):

class TestQueryDistanceTraveled(unittest.TestCase):

class TestQuerySignals(unittest.TestCase):

class TestQuerySpeedLimits(unittest.TestCase):
    
class TestQueryBearings(unittest.TestCase):
    
class TestGetRasterFiles(unittest.TestCase):
    
class TestReprojectRasterfiles(unittest.TestCase):

class TestQueryElevationSeries(unittest.TestCase):
    
class TestSmoothElevation(unittest.TestCase):
    
class TestCombineLidarData(unittest.TestCase):
    
class TestCalculateGrades(unittest.TestCase):
    
class TestInterpolateGeometrySeries(unittest.TestCase):

class TestInterpolateDistanceTraveled(unittest.TestCase):

#--- Route

class TestRouteClass(unittest.TestCase):
    # TestToString
    # TestIntersperseList
    # TestSaveToJSON - saveload?
    # TestQueryPoint
    # TestToGDF

class TestEncodeSeries(unittest.TestCase):

class TestDecodeSeries(unittest.TestCase):
    
class TestEncodeGeometry(unittest.TestCase):
    
class TestDecodeGeometry(unittest.TestCase):
    
class TestLoadFromJson(unittest.TestCase)

if __name__ == '__main__':
    unittest.main()