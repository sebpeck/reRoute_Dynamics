import numpy as np
import shapely

def list_contains(ls, val):
    if (val in ls):
        return True
    else:
        return False
    
def get_bounding_box(shape):
    '''
    get_bounding_box() takes in a shapely shape,
    and then generates a polygon that is a rectangular bounding box.
    
    Parameters:
    shape: any shapely shape
    
    Returns: 
    shapely.polygon that fits around the provided shape
    '''
    # get the bound coordinates
    bbox_cords = shape.bounds
    
    # create a rectangle of those coordinates
    boundary_coords = ((bbox_cords[0], bbox_cords[1]),
                 (bbox_cords[0], bbox_cords[3]),
                 (bbox_cords[2], bbox_cords[3]),
                 (bbox_cords[2], bbox_cords[1]))
    
    # return a shapely polygon of those coordinates
    return shapely.Polygon(boundary_coords)
    

def get_stops(s_gdf, rt_num:str):
    """
    get_stops() takes in an aggregate of bus stops,
    and processes it such that only the stops of the desired route are available
    
    Parameters:
    s_gdf: a geodataframe of bus stop information, including a list of
           routes that arrive there, and geometry of each stop
    rt_num: a string that is the number or ID of the route to be found.
    
    a geodataframe containing solely the data for stops with that route.
    """ 
    
    # Copy the dataframe over so no persistence errors
    temp = s_gdf.copy()
    
    # Fill NAN values with zero
    temp['ROUTE_LIST'] = temp['ROUTE_LIST'].fillna(value = str(0))
    
    # Split the route list into an actual list
    temp['ROUTE_LIST_R'] = temp['ROUTE_LIST'].str.split(" ")
    
    # Use the list_contains method to check if the list contains target route
    temp['RT_MATCH'] = False
    temp['RT_MATCH'] = temp.apply(lambda x: list_contains(ls = list(x['ROUTE_LIST_R']),
                                                          val = rt_num),
                                            axis = 1)
    
    # Filter the stops to only the ones where the target route arrives
    test_stops = temp[(temp['RT_MATCH'] == True) & (temp["IN_SERVICE"] == "Y")]
    return test_stops

def get_lat_lon(pt):
    '''
    get_lat_lon() takes in a shapely geometry point,
    and returns the latitude and longitude as a tuple
    
    Parameters:
    pt: a shapely geometry point
    
    Returns:
    a tuple containing Latitude and Longitude, in that order
    '''
    
    # Conver the point to a string
    
    # Check for formatting
    if ('MULTIPOINT' in str(pt)):
        pt = str(pt)[5:]
        return float((pt)[7:-1].split(" ")[0]), float((pt)[7:-1].split(" ")[1])
    else:
        return float((str(pt)[7:-1]).split(" ")[0]), float((str(pt)[7:-1]).split(" ")[1])


def haversine_formula(lat1:float,
                      lon1:float,
                      lat2:float,
                      lon2:float,
                      r:int = 6371):
    """
    haversine_formula() is a method to calculate the distance
    between two latitude and longitude points on earth, using
    the haversine formula. This has been adapted from:
    
    https://stackoverflow.com/questions/43577086/pandas-calculate-haversine-distance-within-each-group-of-rows
    and
    https://www.omnicalculator.com/other/latitude-longitude-distance
    
    Parameters:
    lat1: latitude of point 1 as a float
    lon1: longitude of point 1 as a float
    lat2: latitude of point 2 as a float
    lon2: longitude of point 2 as a float
    r: radius of earth, default 6371 km

    Returns:
    distance, in kilometers, between the points
    """
    
    # Calculate the first section of the haversine
    a = np.sin((lat2-lat1)/2.0)**2 + \
        np.cos(lat1) * np.cos(lat2) * np.sin((lon2-lon1)/2.0)**2
    
    # Calculate and return the second section of the haversine in km
    return r * 2 * np.arcsin(np.sqrt(a)) / 100 * 1.609344