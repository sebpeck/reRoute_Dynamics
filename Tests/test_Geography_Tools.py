'''
test_Geography_Tools.py

tests the methods in Geography_Tools.py

Notes:
1/14/2025 - I'm uncertain wether or not what I'm doing in TestCompassHeading is 'good coding practice',
            but if I was truly concerned about that, there would be other issues going on in this code.
'''

import unittest
import sys
import os
sys.path.append('../src/reRoute_Dynamics_Core')
import shapely
import pandas as pd
import Geography_Tools as gt
import numpy as np
import geopandas as gpd

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
    

class TestInterpolatePoints(unittest.TestCase):
    # doesnt have handles for negative max_dist. 
    same = (shapely.Point(47.66131, -122.27157), shapely.Point(47.66131, -122.27157))
    within = (shapely.Point(47.66131, -122.27157), shapely.Point(47.66132, -122.27157))
    long_road = (shapely.Point(47.755934, -122.345852), shapely.Point(47.74876, -122.34564)) #797 meters
    def test_same(self):
        result = gt.interpolate_points(*self.same)
        expected = [self.same[0]]
        self.assertEqual(expected, result)
        
    def test_within_bounds(self):
        result = gt.interpolate_points(*self.within)
        expected = [self.within[0]]
        self.assertEqual(expected, result)
    
    def test_one_added(self):
        params = self.long_road + (380,)
        result = len(gt.interpolate_points(*params))
        expected = 2
        
        self.assertEqual(expected, result)
        
    def test_eighthund_added(self):
        params = self.long_road + (1,)
        result = len(gt.interpolate_points(*params))
        expected = 797
        
        self.assertEqual(expected, result)

class TestRepeatIDRemover(unittest.TestCase):
    mt_seq = []
    no_reps = [-1,-1,-1,-1,-1,-1,-1]
    all_reps = [5,5,5,5,5,5,5]
    multi_reps = [2,5,2,5,2,5,3]
    uniques = [5, -1,-1,2,-1,3,-1]
    
    def test_empty_seq(self):
        result = gt.repeat_id_remover(self.mt_seq)
        expected = []
        self.assertEqual(expected,result)
        
    def test_no_reps(self):
        result = gt.repeat_id_remover(self.no_reps)
        expected = self.no_reps
        self.assertEqual(expected,result)
        
    def test_all_reps(self):
        result = gt.repeat_id_remover(self.all_reps)
        expected = [5,-1,-1,-1,-1,-1,-1]
        self.assertEqual(expected, result)
        
    def test_multi_reps(self):
        result = gt.repeat_id_remover(self.multi_reps)
        expected = [2, 5, -1, -1, -1,-1,3]
        self.assertEqual(expected, result)
        
    def test_uniques(self):
        result = gt.repeat_id_remover(self.uniques)
        expected = self.uniques
        self.assertEqual(expected, result)

class TestInterpolateGeometrySeries(unittest.TestCase):
    geo_series = pd.read_csv("../Examples/KC_Example_Data/Example_Route_Data/example_route_table.csv")
    geo_series = pd.Series(geo_series.apply(lambda x: shapely.Point(x.shp_pt_lat, x.shp_pot_lon), axis=1))
    max_dist = 9 #meters
    
    def test_works(self):
        result = gt.interpolate_geometry_series(self.geo_series, self.max_dist)
        self.assertEqual(430, len(result))
        
        

class TestQueryStops(unittest.TestCase):
    geo_series = pd.read_csv("../Examples/KC_Example_Data/Example_Route_Data/example_route_table.csv")
    geo_series = pd.Series(geo_series.apply(lambda x: shapely.Point(x.shp_pt_lat, x.shp_pot_lon), axis=1))
    geo_series = gt.interpolate_geometry_series(geo_series, 20)
    stop_table_path = "../Examples/KC_Example_Data/Example_Route_Data/example_stop_table.csv"

    
    def test_works(self):
        result = gt.query_stops(self.geo_series, self.stop_table_path)
        self.assertEqual(18, len(result) - result.count(-1))
        

class TestQuerySignals(unittest.TestCase):
    # Since this is essentially the same method as querystops, it should be able to work with the same
    # exact test. 
    geo_series = pd.read_csv("../Examples/KC_Example_Data/Example_Route_Data/example_route_table.csv")
    geo_series = pd.Series(geo_series.apply(lambda x: shapely.Point(x.shp_pt_lat, x.shp_pot_lon), axis=1))
    geo_series = gt.interpolate_geometry_series(geo_series, 20)
    stop_table_path = "../Examples/KC_Example_Data/Example_Route_Data/example_stop_table.csv"

    
    def test_works(self):
        result = gt.query_stops(self.geo_series, self.stop_table_path)
        self.assertEqual(18, len(result) - result.count(-1))


class TestQueryBearings(unittest.TestCase):
    geo_series = pd.read_csv("../Examples/KC_Example_Data/Example_Route_Data/example_route_table.csv")
    geo_series = pd.Series(geo_series.apply(lambda x: shapely.Point(x.shp_pt_lat, x.shp_pot_lon), axis=1))
    
    def test_length(self):
        result = len(gt.query_bearings(self.geo_series))
        expected = len(self.geo_series)
        self.assertEqual(expected, result)

        
class TestGetRasterFiles(unittest.TestCase):
    path = "../Examples/KC_Example_Data/Example_Raster_Data/"
    path_with = "../Examples/KC_Example_Data/Example_Raster_Data/Redmond_2014/dtm/"
    def test_no_tif(self):
        result = gt.get_rasterfiles(self.path)
        self.assertEqual(0, len(result))
        
    def test_tif_with_others(self):
        result = gt.get_rasterfiles(self.path_with) # < There should be 4 raster files here
        self.assertEqual(4, len(result))
        

class TestReprojectRasterfiles(unittest.TestCase):
    redmond_rasters_path = "../Examples/KC_Example_Data/Example_Raster_Data/Redmond_2014/dtm/"
    redmond_rasters = gt.get_rasterfiles(redmond_rasters_path)
    
    
    def test_4326(self):
        new_projections = gt.reproject_rasterfiles(self.redmond_rasters)
        self.assertEqual(len(self.redmond_rasters), len(new_projections))
        self.assertEqual(8,len(os.listdir(self.redmond_rasters_path)))
        

class TestQueryElevationSeries(unittest.TestCase):
    redmond_rasters_path = "../Examples/KC_Example_Data/Example_Raster_Data/Redmond_2014/dtm/"
    redmond_rasters = gt.get_rasterfiles(redmond_rasters_path)
    new_projections = gt.reproject_rasterfiles(redmond_rasters)
    pt1 = shapely.Point(47.6778, -122.1574) # close, tract 1
    pt2 = shapely.Point(47.678, -122.1575) # close, tract 1
    pt3 = shapely.Point(47.679, -122.1579) # far, tract 1
    pt4 = shapely.Point(47.6817, -122.1571) # deep, tract 1
    pt5 = shapely.Point(47.6687, -122.1083) # far, onramp, tract 2
    pt6 = shapely.Point(47.6702, -122.1073) # far, bridge of onramp, tract 2
    
    def test_same_tract_close(self):
        pts = [self.pt1, self.pt2, self.pt3]
        result = gt.query_elevation_series(pts, self.new_projections)
        self.assertEqual(3, len(result))
        self.assertEqual(True, result[0]>result[1])
        self.assertEqual(True, result[1]>result[2])
        
    def test_multiple_tracts(self):
        pts = [self.pt1, self.pt4, self.pt6]
        result = gt.query_elevation_series(pts, self.new_projections)
        self.assertEqual(3, len(result))
        self.assertEqual(True, result[0]>result[1])
        self.assertEqual(True, result[0]>result[2])
    

class TestSmoothElevation(unittest.TestCase):
    pts = pd.Series([shapely.Point(47.667,-122.1279),
           shapely.Point(47.6672, -122.1276),
           shapely.Point(47.6673, -122.1255),
           shapely.Point(47.6673, -122.1172),
           shapely.Point(47.6674, -122.111),
           shapely.Point(47.6687,-122.1084),
           shapely.Point(47.6676, -122.1075),
           shapely.Point(47.6702, -122.107)]) # A highway with two bridges at index 1 and index 7

    int_pts = gt.interpolate_geometry_series(pts, 10)
    redmond_rasters_path = "../Examples/KC_Example_Data/Example_Raster_Data/Redmond_2014/dtm/"
    redmond_rasters = gt.get_rasterfiles(redmond_rasters_path)
    new_projections = gt.reproject_rasterfiles(redmond_rasters)
    elevs = gt.query_elevation_series(int_pts, new_projections)
    
    def test_works(self):
        result = gt.smooth_elevation(self.elevs, 20, 3)
        # tTHIS DOES NOT TEST IF THE ELEVATION MAKES REASONABLE SENSE!!!! ONLY WORKS!
        #results = pd.DataFrame([pd.Series(list(self.int_pts)), pd.Series(list(result)), pd.Series(list(self.elevs))]).T
        self.assertEqual(len(self.elevs), len(result))
        self.assertEqual(False, np.isnan(result).any())


class TestCalculateGrades(unittest.TestCase):
    # I am assuming no one is passing negative distances
    dx = [1, 1, 1, 1, 1]
    elevations = [0, .5, -.5, -.5 + 7.5/100, -.5 + 7.5/100+2.5/100]
    
    def test_no_clip(self):
        result = gt.calculate_grades(self.dx, self.elevations, clip=False)
        expected = [50.0, -100.0, 7.5, 2.5, 2.5]
        new = []
        for item in result:
            item = round(item, 1)
            new.append(item)
        self.assertEqual(expected, new)
    
    def test_clip(self):
        result = gt.calculate_grades(self.dx, self.elevations, clip=True)
        expected = [7.5, -7.5, 7.5, 2.5, 2.5]
        new = []
        for item in result:
            item = round(item, 1)
            new.append(item)
        self.assertEqual(expected, new)
        

class TestInterpolateGeometrySeries(unittest.TestCase):
    pts = pd.Series([shapely.Point(47.667,-122.1279),
       shapely.Point(47.6672, -122.1276),
       shapely.Point(47.6673, -122.1255),
       shapely.Point(47.6673, -122.1172),
       shapely.Point(47.6674, -122.111),
       shapely.Point(47.6687,-122.1084),
       shapely.Point(47.6676, -122.1075),
       shapely.Point(47.6702, -122.107)]) # A highway with two bridges at index 1 and index 7

    def test_works(self):
        int_pts = gt.interpolate_geometry_series(self.pts, 10)
        expected = 196
        result = len(int_pts)
        self.assertEqual(expected, result) 

#--- Route

class TestRouteClass(unittest.TestCase):
    pts = pd.Series([shapely.Point(47.667,-122.1279),
       shapely.Point(47.6672, -122.1276),
       shapely.Point(47.6673, -122.1255),
       shapely.Point(47.6673, -122.1172),
       shapely.Point(47.6674, -122.111),
       shapely.Point(47.6687,-122.1084),
       shapely.Point(47.6676, -122.1075),
       shapely.Point(47.6702, -122.107)]) # A highway with two bridges at index 1 and index 7

    int_pts = gt.interpolate_geometry_series(pts, 10)
    redmond_rasters_path = "../Examples/KC_Example_Data/Example_Raster_Data/Redmond_2014/dtm/"
    redmond_rasters = gt.get_rasterfiles(redmond_rasters_path)
    new_projections = gt.reproject_rasterfiles(redmond_rasters)
    elevs = gt.query_elevation_series(int_pts, new_projections)
    
    def test_make_and_methods(self):
        route = gt.Route(self.int_pts, self.elevs)

        self.assertEqual(len(self.elevs), len(route.elevation))
        routesavepath = route.save_to_json('../Examples/KC_Example_Data/Example_Route_Save/redmond_road.json')
        loaded_route = gt.load_from_json(routesavepath)
        self.assertEqual(len(self.elevs),len(loaded_route.elevation))
        pt=route.query_point(0)
        self.assertEqual(self.elevs[0],pt['elevation'])
        self.assertEqual(type(gpd.GeoDataFrame()), type(route.to_gdf()))
        

        
    
'''

    Skipping this one for now since it's not critical to the function.
class TestVerboseLineUpdater(unittest.TestCase):

# I'm gonna be honest, this does not make a whole lot of sense as an independent function. 
class TestQueryElevationChanges(unittest.TestCase):

class TestQueryDistanceTraveled(unittest.TestCase):


This one feels insurmountable right now. How do I create fake geospatial speed limit data?
Do i assemble a bunch of linestrings with different limi... yeah thats the unfortunately tedious solution
hahah frick

class TestQuerySpeedLimits(unittest.TestCase):
    
Coming soon, to a unit test near you...
...EXCEPT YOU WON'T BECAUSE IT's SECRET [read, private methods]!

class TestInterpolateDistanceTraveled(unittest.TestCase):

class TestCombineLidarData(unittest.TestCase):

class TestEncodeSeries(unittest.TestCase):

class TestDecodeSeries(unittest.TestCase):
    
class TestEncodeGeometry(unittest.TestCase):
    
class TestDecodeGeometry(unittest.TestCase):
    
class TestLoadFromJson(unittest.TestCase) < - Already low key tested this one with the Route.
'''
if __name__ == '__main__':
    unittest.main()