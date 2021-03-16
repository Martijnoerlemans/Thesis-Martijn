# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:32:42 2021

@author: Martijn Oerlemans
"""
import math

from h3 import h3
import pandas as pd
import folium
from folium import Map, Marker, GeoJson
import json
from geojson.feature import *
from felyx_gcp_utils.get_gcp_secrets import access_secret_version
from sqlalchemy import create_engine
import os
import psycopg2
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import branca
#Determine granularity of service areas
max_res=15
list_hex_edge_km = []
list_hex_edge_m = []
list_hex_perimeter_km = []
list_hex_perimeter_m = []
list_hex_area_sqkm = []
list_hex_area_sqm = []

for i in range(0,max_res + 1):
    ekm = h3.edge_length(resolution=i, unit='km')
    em = h3.edge_length(resolution=i, unit='m')
    list_hex_edge_km.append(round(ekm,3))
    list_hex_edge_m.append(round(em,3))
    list_hex_perimeter_km.append(round(6 * ekm,3))
    list_hex_perimeter_m.append(round(6 * em,3))
    
    akm = h3.hex_area(resolution=i, unit='km^2')
    am = h3.hex_area(resolution=i, unit='m^2')
    list_hex_area_sqkm.append(round(akm,3))
    list_hex_area_sqm.append(round(am,3))

df_meta = pd.DataFrame({"edge_length_km" : list_hex_edge_km,
                        "perimeter_km" : list_hex_perimeter_km,
                        "area_sqkm": list_hex_area_sqkm,
                        "edge_length_m" : list_hex_edge_m,
                        "perimeter_m" : list_hex_perimeter_m,
                        "area_sqm" : list_hex_area_sqm
                       })
                      
df_meta[["edge_length_km","perimeter_km","area_sqkm", "edge_length_m", "perimeter_m" ,"area_sqm"]]




lat_centr_point = 52.354675
lon_centr_point = 4.875726
list_hex_res = []
list_hex_res_geom = []
list_res = range(0,max_res+1)

for resolution in range(0,max_res + 1):
    #index the point in the H3 hexagon of given index resolution
    h = h3.geo_to_h3(lat=lat_centr_point,lng=lon_centr_point, resolution=resolution)
    list_hex_res.append(h)
    #get the geometry of the hexagon and convert to geojson
    h_geom = { "type" : "Polygon",
               "coordinates": 
                    [h3.h3_to_geo_boundary(h=h,geo_json=True)]
              }
    list_hex_res_geom.append(h_geom)

    
df_resolution_example = pd.DataFrame({"res" : list_res,
                                      "hex_id" : list_hex_res,
                                      "geometry": list_hex_res_geom 
                                     })
df_resolution_example["hex_id_binary"] = df_resolution_example["hex_id"].apply(lambda x: bin(int(x,16))[2:])

pd.set_option('display.max_colwidth',63)
df_resolution_example

map_example = Map(location= [52.354675, 4.875726], zoom_start=5.5, tiles="cartodbpositron", 
                attr= '© <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors © <a href="http://cartodb.com/attributions#basemaps">CartoDB</a>' 
            )
list_features = []
for i,row in df_resolution_example.iterrows():
    feature = Feature(geometry = row["geometry"], id=row["hex_id"], properties = {"resolution": int(row["res"])})
    list_features.append(feature)

feat_collection = FeatureCollection(list_features)
geojson_result = json.dumps(feat_collection)


GeoJson(
        geojson_result,
        style_function=lambda feature: {
            'fillColor': None,
            'color': "green" if feature['properties']['resolution'] %2 == 0 else "red",
            'weight': 2,
            'fillOpacity': 0
        }, 
        name = "Example" 
    ).add_to(map_example)

map_example.save('1_resolutions.html')

#load data from sql database
os.environ["GOOGLE_APPLICATION_CREDENTIALS"]="C:/Users/Martijn Oerlemans/Desktop/felyx-ai-machine-learning-88e8adf81f5b.json"
cnx = create_engine(access_secret_version('felyx-ai-machine-learning','DATABASE_URL_DATAWAREHOUSE', "latest"))


#Data regarding weather, filter out some columns
weather_record_ams = pd.read_sql_query('''SELECT * 
                                       FROM weather_record
                                       WHERE location_id = 1''',cnx)
#data regarding reservations, add columns needed
reservation_ams = pd.read_sql_query('''SELECT *
                                       FROM reservation
                                       WHERE location_id = 1
                                       LIMIT 1000;''',cnx)
        
full_data_ams= pd.read_sql_query('''SELECT a.vehicle_id, a.location_id,
                                  a.reservation_start_time, a.reservation_end_time,
                                  a.start_latitude,a.start_longitude,
                                  a.end_latitude, a.end_longitude,
                                  b.date, b.feels_like_temp, b.snow_level,
                                  b.sun_hours, b.uv_index, b.wind_speed,
                                  b.precipitation, b.humidity, b.visibility,
                                  b.heat_index, b.hour
                                  FROM reservation as a
                                  INNER JOIN weather_record  as b
                                  ON a.location_id = b.location_id
                                  AND a.location_id =1
                                  AND a.location_id_start = 1
                                  AND DATE_TRUNC('day', a.reservation_start_time) = DATE_TRUNC('day', b.date)
                                  AND DATE_PART('hour', a.reservation_start_time) = b.hour
                                  AND a.rent_start_successful
                                  AND NOT a.dev_account
                                  AND (a.rent_end_successful OR a.net_price > 0)
                                  AND a.reservation_end_time < '2021-03-04'
                                  AND a.start_longitude > 4.7''',cnx)

rides = pd.read_sql_query('''SELECT timezone('Europe/Amsterdam', timezone('UTC', a.rent_start_time)) as rent_start_time,
timezone('Europe/Amsterdam', timezone('UTC', a.reservation_start_time)) as reservation_start_time,
timezone('Europe/Amsterdam', timezone('UTC', a.reservation_end_time)) as reservation_end_time,
a.start_latitude,
a.start_longitude,
a.end_latitude,
a.end_longitude,
a.vehicle_id,
a.start_battery_level,
coalesce(sa_start.custom_name,sa_start.default_name) as service_area_start,
b.sun_hours, b.uv_index, b.wind_speed,
b.precipitation, b.humidity, b.visibility,
b.heat_index, b.hour
FROM reservation as a
INNER JOIN weather_record  as b ON a.location_id = b.location_id
LEFT JOIN service_area sa_start on st_contains(sa_start.wgs84_polygon,ST_SetSRID(st_makepoint(a.start_longitude, a.start_latitude),4326)) AND sa_start.activated
LEFT JOIN service_area sa_end on st_contains(sa_end.wgs84_polygon,ST_SetSRID(st_makepoint(a.end_longitude, a.end_latitude),4326)) AND sa_end.activated
WHERE a.location_id = 1
AND DATE_TRUNC('day', reservation_start_time) = DATE_TRUNC('day', b.date)
AND DATE_PART('hour', reservation_start_time) = b.hour
and reservation_start_time < '2021-03-11' 
AND rent_start_time is not null
AND NOT a.dev_account
AND a.rent_start_successful
AND sa_start.id IS NOT NULL
AND sa_end.id IS NOT NULL
AND a.start_longitude > 4.7''', cnx)

full_data_ams = rides
trainstations_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/trainstations_dwh.xlsx')
trainstations_dwh_ams = trainstations_dwh[trainstations_dwh['location_id']==1]
uni_hbo_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/uni_hbo_dwh.xlsx')
uni_hbo_dwh_ams = uni_hbo_dwh[uni_hbo_dwh['location_id']==1]

#full_data_ams = full_data_ams[full_data_ams['location_id_start']==1]
#resolution wanted
service_area_resolution = 9
#converting long,lat columns to hexagons
full_data_ams['start_hexagon'] = full_data_ams.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['start_latitude'],
                                    lng=row['start_longitude'], 
                                    resolution=service_area_resolution), axis=1)

full_data_ams['end_hexagon'] = full_data_ams.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['end_latitude'],
                                    lng=row['end_longitude'], 
                                    resolution=service_area_resolution), axis=1)

full_data_ams_100= full_data_ams.head(100)
trainstations_dwh_ams['hexagon'] = trainstations_dwh_ams.apply(lambda row: 
                            h3.geo_to_h3(lat=row['latitude'],
                            lng=row['longitude'], 
                            resolution=service_area_resolution), axis=1)
uni_hbo_dwh_ams['hexagon'] = uni_hbo_dwh_ams.apply(lambda row: 
                            h3.geo_to_h3(lat=row['latitude'],
                            lng=row['longitude'], 
                            resolution=service_area_resolution), axis=1)
#save data to excel file
#full_data_ams.to_excel("full_data_ams.xlsx")
#1832138
data = full_data_ams.groupby(full_data_ams.reservation_start_time.hour)

data =full_data_ams.groupby([pd.Grouper(key='reservation_start_time',freq='H'),
                             'start_hexagon']).agg({'sun_hours' : ['first'],
                            'reservation_start_time' : ['count'],'service_area_start':'first',
                            'start_battery_level' : ['mean'],'uv_index':['first'],
                            'wind_speed':['first'],'precipitation':['first'],
                            'humidity':['first'],'visibility':['first'],'heat_index':['first']
                            ,}).reset_index()
data100 = data.head(100)
unique_hexagons_ams_start = full_data_ams.start_hexagon.unique()
unique_hexagons_ams_end = full_data_ams.end_hexagon.unique()

random_hexagon ='89196953157ffff'
# function to get boundaries of a hexagon in latitudes and longitudes
h3.h3_to_geo_boundary(random_hexagon)

def visualize_hexagons(hexagons, color, folium_map=None):
    """
    hexagons is a list of hexcluster. Each hexcluster is a list of hexagons. 
    eg. [[hex1, hex2], [hex3, hex4]]
    """
    polylines = []
    lat = []
    lng = []
    for hex in hexagons:
        polygons = h3.h3_set_to_multi_polygon([hex], geo_json=False)
        # flatten polygons into loops.
        outlines = [loop for polygon in polygons for loop in polygon]
        polyline = [outline + [outline[0]] for outline in outlines][0]
        lat.extend(map(lambda v:v[0],polyline))
        lng.extend(map(lambda v:v[1],polyline))
        polylines.append(polyline)
    
    if folium_map is None:
        m = folium.Map(location=[sum(lat)/len(lat), sum(lng)/len(lng)], zoom_start=13, tiles='cartodbpositron')
    else:
        m = folium_map
    for polyline in polylines:
        my_PolyLine=folium.PolyLine(locations=polyline,weight=0.2,color=color,fill=True,fill_color = color,fill_opacity = 0.2)
        m.add_child(my_PolyLine)
    return m

m = visualize_hexagons(unique_hexagons_ams_start,
                               'green' 
                               ,None)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\test.html")



def colorFader(c1,c2,mix): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

#calculates haversine distances between two long lat points in meters

def haversine(coord1, coord2):
    R = 6372800  # Earth radius in meters
    lat1, lon1 = coord1
    lat2, lon2 = coord2
    
    phi1, phi2 = math.radians(lat1), math.radians(lat2) 
    dphi       = math.radians(lat2 - lat1)
    dlambda    = math.radians(lon2 - lon1)
    
    a = math.sin(dphi/2)**2 + \
        math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
    
    return 2*R*math.atan2(math.sqrt(a), math.sqrt(1 - a))

def minimum_distance_to_multiple_points(latitudes,longitudes,latitude,longitude):
    for i in range(len(latitudes)):
        minimum_distance = 100000000
        distance_i = haversine([latitudes[i],longitudes[i]],[latitude,longitude])
        if  distance_i < minimum_distance:
            minimum_distance = distance_i
            minimum_distance_coordinates = [latitudes[i],longitudes[i]]
    return minimum_distance, minimum_distance_coordinates

full_data_ams['distance_closest_trainstation'] = full_data_ams.apply(lambda row: 
                                    minimum_distance_to_multiple_points(
                                        trainstations_dwh_ams['latitude'],trainstations_dwh_ams['longitude']
                                        ,row['start_latitude'],row['start_longitude'])[0]
                                    , axis=1)
    
full_data_ams_scooters_used =full_data_ams['start_hexagon'].value_counts()
maximum_scooters_used = max(full_data_ams_scooters_used)


#all reservations in Amsterdam
c1 = 'green'
c2='red'

for i in range(len(full_data_ams_scooters_used.unique())):
    if i == 0:
        m = visualize_hexagons(full_data_ams_scooters_used.index[full_data_ams_scooters_used 
                               == full_data_ams_scooters_used.unique()[0]],
                               colorFader(c1,c2,
                                full_data_ams_scooters_used.unique()[0]/maximum_scooters_used) 
                               ,None)
    else:
        m = visualize_hexagons(full_data_ams_scooters_used.index[full_data_ams_scooters_used
                               == full_data_ams_scooters_used.unique()[i]],
                               colorFader(c1,c2,
                                full_data_ams_scooters_used.unique()[i]/maximum_scooters_used) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_scooters_used)

colormap = colormap.to_step(index=[0, maximum_scooters_used* (1/10), maximum_scooters_used* (2/10)
                                   , maximum_scooters_used* (3/10), maximum_scooters_used* (4/10)
                                   , maximum_scooters_used* (5/10), maximum_scooters_used* (6/10)
                                   , maximum_scooters_used* (7/10), maximum_scooters_used* (8/10),
                                   maximum_scooters_used* (9/10), maximum_scooters_used])
colormap.caption = 'Scooters reserved and used per hexagon in Amsterdam (service_area)'
colormap.add_to(m)
folium.Marker(location=[52.920731882528166, 4.670521006248073],popup='Total scooters reserved in Amsterdam :' 
              + str(len(full_data_ams)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\scooters_used_amsterdam.html")