# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:32:42 2021

@author: Martijn Oerlemans
"""

from h3 import h3
import pandas as pd
from folium import Map, Marker, GeoJson
import json
from geojson.feature import *
from felyx_gcp_utils.get_gcp_secrets import access_secret_version
from sqlalchemy import create_engine
import os
import psycopg2
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
        
uni_hbo_dwh = pd.read_sql_query('''SELECT *
                                FROM uni_hbo_dwh''',cnx)
uni_hbo_dwh_ams = pd.read_sql_query('''SELECT *
                                       FROM uni_hbo_dwh
                                       WHERE location_id = 1
                                       LIMIT 1000;''',cnx)
treinstations_dwh_ams = pd.read_sql_query('''SELECT *
                                       FROM trainstations_dwh
                                       WHERE location_id = 1
                                       LIMIT 1000;''',cnx)

full_data_ams = pd.read_sql_query('''SELECT a.*, b.*
                                  FROM reservation as a
                                  INNER JOIN weather_record  as b
                                  ON a.location_id = b.location_id
                                  AND a.location_id =1 -- Amsterdam
                                  AND DATE_TRUNC('day', a.reservation_start_time) = DATE_TRUNC('day', b.date)
                                  AND DATE_PART('hour', a.reservation_start_time) = b.hour
                                  AND a.rent_start_successful
                                  AND NOT a.dev_account
                                  AND (a.rent_end_successful OR a.net_price > 0)
                                  AND a.reservation_end_time < '2021-03-04' ''',cnx)
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
#save data to excel file
#full_data_ams.to_excel("full_data_ams.xlsx")
#1832138
unique_hexagons_ams_start = full_data_ams.start_hexagon.unique()
unique_hexagons_ams_end = full_data_ams.end_hexagon.unique()