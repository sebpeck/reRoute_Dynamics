'''
test_Geography_Tools.py

tests the methods in Geography_Tools.py

Notes:
1/14/2025 - I'm uncertain wether or not what I'm doing in TestCompassHeading is 'good coding practice',
            but if I was truly concerned about that, there would be other issues going on in this code.
'''

import unittest
import sys
sys.path.append('../src/reRoute_Dynamics_Core')
import shapely
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
    def test_north(self):
        self.assertEqual('N', gt.compass_heading(0))
        self.assertEqual('N', gt.compass_heading(360))
        self.assertEqual('N', gt.compass_heading(720))
        
    def test_non_standard(self):
        self.assertEqual('SW', gt.compass_heading(180+42))
    
    def test_pass_back_and_forth(self):
        self.assertEqual('E', gt.compass_heading(gt.heading_to_angle('E')))

class TestHeadingToAngle(unittest.TestCase):
    heads = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
    expects = [0, 45, 90, 90+45, 180, 180+45, 180+90, 180+90+45]
    def test_all(self):
        for i in range(len(self.heads)-1):
            self.assertEqual(self.expects[i], gt.heading_to_angle(self.heads[i]))

class TestPointBearing(unittest.TestCase):
    
    same = (47.66131, -122.27157, 47.66131, -122.27157)
    short_north = (47.66131, -122.27157, 47.66305, -122.27139, "Compass")
    long_north = (43.03527, -122.27157, 47.66305, -122.27139, "Compass")
    def test_same_point(self):
        result = gt.point_bearing(*self.same)
        expected = 0
        self.assertEqual(expected, result)
        
    def test_short_north(self):
        result = gt.point_bearing(*self.short_north)
        expected = 'N'
        self.assertEqual(expected, result)
    
    def test_long_north(self):
        result = gt.point_bearing(*self.long_north)
        expected = 'N'
        self.assertEqual(expected, result)

class TestBoundingBox(unittest.TestCase):
    point = shapely.Point(0, 0)
    line = shapely.LineString([(-1, 0), (1, 0)])
    triangle = shapely.Polygon([(-1,0),(0,1),(1,0), (-1,0)])
    goofy = shapely.LineString([(-1, 0),(5, 2),(-3, 7.5),(8, 5),(7, 10)])
    
    def test_point(self):
        result = gt.get_bounding_box(self.point)
        expected = shapely.Polygon(((0,0), (0,0), (0,0), (0,0)))
        self.assertEqual(expected, result)
        
    def test_line(self):
        result = gt.get_bounding_box(self.line)
        expected = shapely.Polygon(((-1,0), (-1,0), (1,0), (1,0)))
        self.assertEqual(expected, result)
        
    def test_triangle(self):
        result = gt.get_bounding_box(self.triangle)
        expected = shapely.Polygon(((-1,0), (-1,1), (1,1), (1,0)))
        self.assertEqual(expected, result)
        
    # Super, Hexagon.
    
    def test_goofy(self):
        result = gt.get_bounding_box(self.goofy)
        expected = shapely.Polygon(((-3,0), (-3,10), (8,10), (8,0)))
        self.assertEqual(expected, result)
    
'''
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
'''
if __name__ == '__main__':
    unittest.main()