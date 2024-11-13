'''
Geography_Tools.py
'''

def haversine_formula(x1, y1, x2, y2):
    '''
    haversine_formula takex the latitude and longitude of two separate points,
    and calculates the distance in kilometers between those two points as a crow flies.
    
    Params:
    x1 - latitude pt 1 in degrees
    y1 - longitude pt 1 in degrees
    x2 - latitude pt 2 in degrees
    y2 - longitude pt 2 in degrees
    
    Returns:
    distance between the points in kilometers.
    '''
    
    # convert the degree coordinates to radians
    lat1=np.radians(x1)
    lon1=np.radians(y1)
    lat2=np.radians(x2)
    lon2=np.radians(y2)
    
    # Get the radius of earth in kilometers
    R = 6373.0
    
    # calculate the difference in radians
    dlon=lon2-lon1
    dlat=lat2-lat1
    
    # calculate the distance using the haversine formula.
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))

    distance = R * c
    
    
    # return the distance in km.
    return distance


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

### ---- SERIES QUERY METHODS ---- ###

def query_stops(geometry, stop_table_path, key='stop_id', epsg_from = 4326, epsg_to=4326, margin=10):
    '''
    query_stops takes in a geometry point series for a given route, and a path to
    all stop ids and geometry, and then returns a list of stop id's in the order the
    bus will arrive.
    
    Params:
    geometry - a series of shapely points representing longitude and latitude
    stop_table_path - a path to the dataset containing stop geodata csv, as str
    key - column identifier for the stop id (or whatever value you want) as str
    epsg_from - the epsg projection the csv is in as an int, default is 4326
    epsg_to - the epsg projection the data should be reporojected to as an int, default is 4326
    margin - distance, in meters, that the stop is allowed to be from a given point to qualify, as int.
    
    Returns:
    list of stop ids (or other queried values from the data), with invalid points being labeled -1.
    '''
    
    # Check the shape boundaries for the geometric series
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))
    
    # read the data.
    stop_table = pd.read_csv(stop_table_path)
    
    # convert the data to a shapely geometry pandas series
    stop_table['geometry'] = stop_table.apply(lambda x: shapely.Point(x.stop_lon, x.stop_lat), axis=1)
    
    # convert the data to a geodataframe, and re-set the crs to the specified
    stop_table = gpd.GeoDataFrame(stop_table).set_crs(epsg = epsg_from).to_crs(epsg = epsg_to)
    
    # filter the stop table to only be the stops within the boundary
    stop_table = stop_table[stop_table['geometry'].apply(lambda x: shape_bounds.contains(x)) == True]
    
    
    # check for stops within a set distance of each geometry point
    stop_id_list = []
    
    # loop through the geimetry
    for point in list(geometry):
        
        # re-set the pont name because i dont want to re-change them all
        route_pt = point
        
        # get the distances using the haversine formula, converted to meters
        stop_table['dist'] = stop_table.geometry.apply(lambda x: haversine_formula(x.y, x.x, route_pt.y, route_pt.x)*1000)
        
        # filter the data to yeild the closest stop
        closest_stop = stop_table[stop_table['dist'] == stop_table['dist'].min()].reset_index().iloc[0]
        
        # if the closest stop is within the margin, append it to the list. Otherwise, append -1.
        stop_id_list.append(int(closest_stop['dist']<margin)*closest_stop[key] - int(closest_stop['dist']>margin))
        
    # return the stop id list.
    return stop_id_list


def query_distance_traveled(geometry_series):
    '''
    query_distance_traveled takes in a series of points,
    and then calculates the distance traveled from point to point, in kilometers.

    Parameters:
    geometry_series - a series of shapely points of lat and lon.

    Returns:
    an iterable representing the change in distance between each point and the next.
    '''

    # get the current geometry series
    currents = geometry_series

    # shift the geometry series up by one, which will be the next ponts.
    nexts = geometry_series.shift(1)

    # combine the two into a dataframe.
    data = pd.concat([currents, nexts], axis=1, keys=['current', 'next']).reset_index(drop=True)

    # get the initial point and make sure it's the first value.
    initial_point = list(data['current'])[0]

    # set the nan value in the beginning of the 'nexts' to the initial, representing no distance traveled.
    data.loc[0, 'next'] = initial_point

    # apply the haversine formula to the data, which will yeild a series of changes in distance.
    dXs = data.apply(lambda x: haversine_formula(x.current.y, x.current.x, x.next.y, x.next.x), axis=1)

    # return the series of the changes in distance, in km
    return dXs.values


def query_signals(geometry, signal_data_path, key="SIGNAL_ID", epsg = 4326, margin = 10):
    '''
    query_signals takes in a geometric series of shapely points,
    a path to stoplight signal data, and returns a list of the signal id at
    each point in that series.
    
    Parameters:
    geometry - a series of shapely points of lat and lon.
    signal_data_path - path, as str, to a shapefile containing stop signal data.
    key - the column identifier of signal id (or other target value) to be queried
    epsg - the epsg the data will be re-set as. Default is 4326.
    margin - the distance, in meters, that the signal is allowed to be from a given point to qualify.
    
    Returns:
    list of signal IDs, or other info about the signal, as specified.
    -1 if none found.
    '''
    
    # get the bounding box of the geometry
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))
    
    # Get the signal data.
    signals = gpd.read_file(signal_data_path).to_crs(epsg = epsg)
    
    # get the signals that are within the bounds of the data.
    signals = signals[signals.apply(lambda x: shape_bounds.contains(x.geometry), axis=1).reset_index(drop=True) == True]
    
    # Generate empty list for signal ids
    signal_id_list = []
    
    # loop through the geometry
    for point in list(geometry):
        
        # use the haversine formula to get the distances between the current point and each point in the series, in meters
        signals['dist'] = signals.geometry.apply(lambda x: haversine_formula(x.y, x.x, point.y, point.x)*1000)
        
        # get the closest signal through the minimum distance
        closest_signal = signals[signals['dist'] == signals['dist'].min()].reset_index().iloc[0]
        
        # check if the signal is within the margin, and if so, append it to the list. If not, append -1.
        signal_id_list.append((int(closest_signal['dist'])<margin)*closest_signal[key] - int(closest_signal['dist']>margin))
        
    # return the list.
    return signal_id_list



def query_speed_limits(geometry, limit_data_path, key="SPEED_LIM", margin = 1):
    '''
    query_speed_limits takes a geometric iterable of shapely points,
    a path to speed limit geodata, and returns the corresponding speed limits at
    each point.
    
    Params:
    geometry - a series of shapely points of lat and lon
    limit_data_path - path to shapefile of speed limit data as str
    key - column identifier of a speed limit id or other value as str
    margin - distance a point can be from the street to qualify in meters as int
    
    Returns:
    list of speed limits corresponding to the geometry.
    
    Notes:
    11/13/2024 - This can be adjusted such that instead of the first value, it's the true closest.
                 Also can probably be combined with the other queries in such a way that it reduces individual
                 methods. 
    '''
    # DEPENDING ON THE FILE, THIS CAN TAKE A WHILE.
    
    # get the bounding box of the route.
    shape_bounds = get_bounding_box(shapely.LineString(list(geometry)))
    
    # time benchmark because this takes a while
    print("TP1: {}".format(dt.datetime.now()))
    
    # get the speed limit data from the shapefile
    limits = gpd.read_file(limit_data_path)
    
    # timepoint
    print("TP2: {}".format(dt.datetime.now()))
    
    # get the streets that are within the bounding box
    bound_streets = limits[limits.apply(lambda x: shape_bounds.contains(x.geometry), axis=1).reset_index(drop=True) == True][['geometry', key]]

    # timepoint
    print("TP3: {}".format(dt.datetime.now()))
    
    # set up list.
    limit_list = []
    
    # timepoint
    print("TP4: {}".format(dt.datetime.now()))
    
    # loop through the geometry:
    for point in list(geometry):
        
        # get the speed limit at the point through checking that the distance is within the margin. Use first value that appears.
        limit = bound_streets[bound_streets.geometry.apply(lambda x: x.distance(point)*1000 < margin) == True].reset_index().iloc[0][key]
        
        # append the limit to the list.
        limit_list.append(limit)
    
    # timepoint
    print("TP5: {}".format(dt.datetime.now()))
    
    # while there is a speed limit of zero, swap them out for the next nearest speed limit that isn't zero.
    while (0 in limit_list):
        for i in range(len(limit_list)):
            print(i, end='\r')
            if limit_list[i] == 0:
                try:
                    limit_list[i] = limit_list[i+1]
                except:
                    limit_list[i] = limit_list[i-1]
                    
    # timepoint
    print("TP6: {}".format(dt.datetime.now()))

    # return the limit list.
    return limit_list


### ---- ELEVATION TOOLS ---- ###

