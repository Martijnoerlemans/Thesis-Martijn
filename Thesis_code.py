# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 16:32:42 2021

@author: Martijn Oerlemans
"""
# %% packages
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
import shapefile
import geopandas as gpd
import contextily as ctx
from shapely import wkt
from tqdm import tqdm
from datetime import datetime
import rtree
from shapely.geometry import Polygon
from shapely.geometry import Point
from felyx_gcp_utils.read_write_gcp_storage import write_pickle_to_bucket, load_pickle_from_bucket, load_csv_from_bucket

# %%Code to determine which hexagon size we want to take (subjective)
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
#table with different levels of hexagons
df_meta = pd.DataFrame({"edge_length_km" : list_hex_edge_km,
                        "perimeter_km" : list_hex_perimeter_km,
                        "area_sqkm": list_hex_area_sqkm,
                        "edge_length_m" : list_hex_edge_m,
                        "perimeter_m" : list_hex_perimeter_m,
                        "area_sqm" : list_hex_area_sqm
                       })
                      
df_meta[["edge_length_km","perimeter_km","area_sqkm", "edge_length_m", "perimeter_m" ,"area_sqm"]]



#Visualizing a point in different levels of hexagons
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

# %%loading data from sql database and other sources
#setting up the connection
os.environ["GOOGLE_APPLICATION_CREDENTIALS"]="C:/Users/Martijn Oerlemans/Desktop/felyx-ai-machine-learning-88e8adf81f5b.json"
cnx = create_engine(access_secret_version('felyx-ai-machine-learning','DATABASE_URL_DATAWAREHOUSE', "latest"))


# #Data regarding weather, filter out some columns
# weather_record_ams = pd.read_sql_query('''SELECT * 
#                                        FROM weather_record
#                                        WHERE location_id = 1''',cnx)
# #data regarding reservations, add columns needed
# reservation_ams = pd.read_sql_query('''SELECT *
#                                        FROM reservation
#                                        WHERE location_id = 1
#                                        LIMIT 1000;''',cnx)
        
# full_data_ams= pd.read_sql_query('''SELECT a.vehicle_id, a.location_id,
#                                   a.reservation_start_time, a.reservation_end_time,
#                                   a.start_latitude,a.start_longitude,
#                                   a.end_latitude, a.end_longitude,
#                                   b.date, b.feels_like_temp, b.snow_level,
#                                   b.sun_hours, b.uv_index, b.wind_speed,
#                                   b.precipitation, b.humidity, b.visibility,
#                                   b.heat_index, b.hour
#                                   FROM reservation as a
#                                   INNER JOIN weather_record  as b
#                                   ON a.location_id = b.location_id
#                                   AND a.location_id = 1
#                                   AND a.location_id_start = 1
#                                   AND DATE_TRUNC('day', a.reservation_start_time) = DATE_TRUNC('day', b.date)
#                                   AND DATE_PART('hour', a.reservation_start_time) = b.hour
#                                   AND a.rent_start_successful
#                                   AND NOT a.dev_account
#                                   AND (a.rent_end_successful OR a.net_price > 0)
#                                   AND a.reservation_end_time < '2021-03-04'
#                                   AND a.start_longitude > 4.7''',cnx)

full_data_ams = pd.read_sql_query('''SELECT timezone('Europe/Amsterdam', timezone('UTC', a.rent_start_time)) as rent_start_time,
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
AND a.start_longitude > 4.7
LIMIT 1000;
''', cnx)



trainstations_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/trainstations_dwh.xlsx')
trainstations_dwh_ams = trainstations_dwh[trainstations_dwh['location_id']==1].reset_index()

uni_hbo_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/uni_hbo_dwh.xlsx')
uni_hbo_dwh_ams = uni_hbo_dwh[uni_hbo_dwh['location_id']==1].reset_index()

# bodemgebruik_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/bodemgebruik_dwh.xlsx')
# bodemgebruik_dwh_ams = bodemgebruik_dwh[bodemgebruik_dwh['location_id']==1]

# cbs_zipcode_2019_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/cbs_zipcode_2019_dwh.xlsx')

metro_dwh = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/metro_dwh.xlsx')
metro_dwh_ams = metro_dwh[metro_dwh['location_id']==1].reset_index()

# kwb_2020 = pd.read_excel(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/kwb-2020.xls')
# kwb_2020_ams = kwb_2020[kwb_2020['gm_naam']=='Amsterdam']

#loading in environmental data
#52.290819,4.786194,52.410362,5.038948
#http://bboxfinder.com/#52.290819,4.786194,52.410362,5.038948
# EPSG: 28992 Amersfoort / RD New
bbox = (113995.7727,478262.1749,131303.8192,491450.0742)

fp = r'C:\Users\Martijn Oerlemans\Downloads\2021-cbs_vk500_2020_v1\CBS_vk500_2020_v1.shp'
CBS_vk500_2020_v1 = gpd.read_file(fp,bbox=bbox).to_crs(epsg=4326)

fp = r'C:\Users\Martijn Oerlemans\Downloads\2021-cbs_vk500_2019_v2\CBS_vk500_2019_v2.shp'
CBS_vk500_2019_v2 = gpd.read_file(fp,bbox=bbox).to_crs(epsg=4326)

fp = r'C:\Users\Martijn Oerlemans\Downloads\2021-cbs_vk500_2018_v3\CBS_vk500_2018_v3.shp'
CBS_vk500_2018_v3 = gpd.read_file(fp,bbox=bbox).to_crs(epsg=4326)

fp = r'C:\Users\Martijn Oerlemans\Downloads\WijkBuurtkaart_2020_v1\gemeente_2020_v1.shp'
gemeente_2020_v1 = gpd.read_file(fp,bbox=bbox).to_crs(epsg=4326)

fp = r'C:\Users\Martijn Oerlemans\Downloads\bestandbodemgebruik2015\BBG2015.shp'
BBG2015 = gpd.read_file(fp,bbox=bbox).to_crs(epsg=4326) 
BBG2015 = BBG2015.drop(columns = ['Shape_Leng', 'Shape_Area'], axis=1)
BBG2015.Omschrijvi.unique()
OmschrijvingToKeep = ['Vliegveld',
       'Detailhandel en horeca', 'Openbare voorziening',
       'Sociaal-culturele voorziening','Park en plantsoen',
       'Sportterrein', 'Volkstuin', 'Dagrecreatief terrein', 
       'Bos', 'Water met recreatieve functie']
BBG2015  = BBG2015[BBG2015['Omschrijvi'].isin(OmschrijvingToKeep)]

BBG2015_head = BBG2015.head()
#covid_data
covid = load_csv_from_bucket('felyx-ai-machine-learning', 'general/', 'COVID_data/', 'covid_variables.csv',";")

#Data regarding operationality
operationality = pd.read_sql_query('''SELECT * 
                                        FROM operationality
                                        WHERE vehicle_State_id = 16
                                        AND location_id = 1''',cnx)
                                        

# %% creating data and appending the dataset
#full_data_ams = full_data_ams[full_data_ams['location_id_start']==1]
#resolution wanted
service_area_resolution = 9
#converting long,lat columns to hexagons

full_data_ams['start_hexagon'] = full_data_ams.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['start_latitude'],
                                    lng=row['start_longitude'], 
                                    resolution=service_area_resolution), axis=1)

    
# full_data_ams['end_hexagon'] = full_data_ams.apply(lambda row: 
#                                     h3.geo_to_h3(lat=row['end_latitude'],
#                                     lng=row['end_longitude'], 
#                                     resolution=service_area_resolution), axis=1)

# full_data_ams_100= full_data_ams.head(100)
# trainstations_dwh_ams['hexagon'] = trainstations_dwh_ams.apply(lambda row: 
#                             h3.geo_to_h3(lat=row['latitude'],
#                             lng=row['longitude'], 
#                             resolution=service_area_resolution), axis=1)
# uni_hbo_dwh_ams['hexagon'] = uni_hbo_dwh_ams.apply(lambda row: 
#                             h3.geo_to_h3(lat=row['latitude'],
#                             lng=row['longitude'], 
#                             resolution=service_area_resolution), axis=1)

    
#converting hexagons to polygons or points
full_data_ams['Polygon']  = full_data_ams.apply(lambda row: 
                                Polygon(h3.h3_to_geo_boundary(row['start_hexagon'],
                                geo_json=True)).wkt, axis=1).apply(wkt.loads)
# full_data_ams['centroid_start_hexagon']=full_data_ams.apply(lambda row: 
#                                 Point(h3.h3_to_geo(row['start_hexagon'])).wkt, axis=1).apply(wkt.loads)    
#set geometry
full_data_ams = full_data_ams.set_geometry(full_data_ams['Polygon'])
  
# full_data_ams = full_data_ams.set_geometry(full_data_ams['centroid_start_hexagon'])
#converting series to geoseries
full_data_ams['Polygon'] = gpd.GeoSeries(full_data_ams['Polygon'],crs = 'EPSG:4326')

# full_data_ams['geometry'] = gpd.GeoSeries(full_data_ams['centroid_start_hexagon'],crs = 'EPSG:4326')
#save data to excel file
#full_data_ams.to_excel("full_data_ams.xlsx")
#1832138
#check what kind of crs
# full_data_ams.crs
# CBS_vk500_2020_v1.crs
#convert to geodataframe
full_data_ams_gdf = gpd.GeoDataFrame(full_data_ams, crs="EPSG:4326",geometry = 'Polygon')
# full_data_ams_gdf = gpd.GeoDataFrame(full_data_ams, crs="EPSG:4326",geometry = 'centroid_start_hexagon')
#geospatial join on hexagons that are in certain polygon of bodemgebruik
full_data_ams_gdf_bbg= gpd.sjoin(full_data_ams_gdf,
                         BBG2015,
                         how="left", 
                         op="overlaps")
#drop irrelevant columns
full_data_ams_gdf_bbg= full_data_ams_gdf_bbg.drop(['index_right','BG2015', 'Hoofdgroep'], 1)
#drop duplicates
full_data_ams_gdf_bbg = full_data_ams_gdf_bbg.drop_duplicates()
#Cluster bodemgebruik per instance
full_data_ams_gdf_groundgroups = full_data_ams_gdf_bbg['Omschrijvi'].groupby([full_data_ams_gdf_bbg.index]).apply(list)
#make into a list
full_data_ams_gdf_groundgroups = full_data_ams_gdf_groundgroups.tolist() # list of lists is required for the following code
#make into dataframe
df = pd.DataFrame(
    {'groups':
full_data_ams_gdf_groundgroups
    }, columns=['groups'])
#create dummy for all bodemgebruik types and handle the NaNs
df =df.groups.apply(lambda x: pd.Series([1] * len(x), index=x)).fillna(0, downcast='infer')
#Concatenate the dummies table to the original table
full_data_ams_bbg = pd.concat([full_data_ams_gdf, df], axis=1)
data = full_data_ams_bbg
#delete this big boy table as it is very large
del BBG2015



full_data_ams_bbg_CBS = gpd.sjoin(full_data_ams_bbg,
                         CBS_vk500_2020_v1,
                         how="left", 
                         op="overlaps")
#convert back to dataframe
full_data_ams_bbg_CBS =pd.DataFrame(full_data_ams_bbg_CBS)
def replace(listofnumbers):
    
    a = [0 if x<1 else x for x in listofnumbers]
    return a

columnsused = CBS_vk500_2020_v1.columns
for i in range(len(columnsused)):
    column_i = columnsused[i]
    
    a = full_data_ams_bbg_CBS[column_i].groupby([full_data_ams_bbg_CBS.index]).apply(list)
    
    if type(a[0][0]) == float:
        a = a.apply(lambda x: replace(x))
        a=a.fillna(0, downcast='infer')
        df = pd.DataFrame(
            {column_i: a}, columns=[column_i])
        df = df[column_i].apply(lambda x: np.mean(x))
        data = pd.concat([data, df], axis=1)


tussenstopdata = data
#loop to create text for function
text = ''
for i in range(len(data.columns)):
    if type(data.columns[i]) == str:
        text = text + ''', \'{}'''.format(data.columns[i]) + '''\' : ['first']'''

#group the data by hexagon and reservation start time per 3 hours
data =data.groupby([pd.Grouper(key='reservation_start_time',freq='3H'),
                             'start_hexagon']).agg({'reservation_start_time' : ['count'],
                                                    'vehicle_id' : ['first'], 'start_battery_level' : ['mean'], 
                                                    'service_area_start' : ['first'], 'sun_hours' : ['first'], 'uv_index' : ['first'], 
                                                    'wind_speed' : ['first'], 'precipitation' : ['first'], 'humidity' : ['first'], 
                                                    'visibility' : ['first'], 'heat_index' : ['first'], 'hour' : ['first'], 
                                                    'Polygon' : ['first'], 'geometry' : ['first'], 
                                                    'Bos' : ['first'], 'Dagrecreatief terrein' : ['first'], 
                                                    'Detailhandel en horeca' : ['first'], 'Openbare voorziening' : ['first'], 
                                                    'Park en plantsoen' : ['first'], 'Sociaal-culturele voorziening' : ['first'], 
                                                    'Sportterrein' : ['first'], 'Volkstuin' : ['first'], 
                                                    'Water met recreatieve functie' : ['first'], 'INWONER' : ['first'], 
                                                    'MAN' : ['first'], 'VROUW' : ['first'], 'INW_014' : ['first'], 
                                                    'INW_1524' : ['first'], 'INW_2544' : ['first'], 'INW_4564' : ['first'], 
                                                    'INW_65PL' : ['first'], 'P_NL_ACHTG' : ['first'], 'P_WE_MIG_A' : ['first'], 
                                                    'P_NW_MIG_A' : ['first'], 'AANTAL_HH' : ['first'], 'TOTHH_EENP' : ['first'], 
                                                    'TOTHH_MPZK' : ['first'], 'HH_EENOUD' : ['first'], 'HH_TWEEOUD' : ['first'], 
                                                    'GEM_HH_GR' : ['first'], 'WONING' : ['first'], 'WONVOOR45' : ['first'], 
                                                    'WON_4564' : ['first'], 'WON_6574' : ['first'], 'WON_7584' : ['first'], 
                                                    'WON_8594' : ['first'], 'WON_9504' : ['first'], 'WON_0514' : ['first'], 
                                                    'WON_1524' : ['first'], 'WON_MRGEZ' : ['first'], 'P_KOOPWON' : ['first'], 
                                                    'P_HUURWON' : ['first'], 'WON_HCORP' : ['first'], 'WON_NBEW' : ['first'], 
                                                    'WOZWONING' : ['first'], 'UITKMINAOW' : ['first'], 'OAD' : ['first'], 
                                                    'STED' : ['first']}).reset_index()
#get rid of the multi index                                                    
data.columns = [col[0] for col in data.columns.values]
#function to count the amount of surrounding hexagons
def Count_Surrounding_Hexagons(hexagon,list_of_hexagons):
    hexagons_used = list(list_of_hexagons)
    surrounding_hexagons =list(h3.k_ring(hexagon,1))[1:]
    count = 0
    for i in range(0,6):
        if surrounding_hexagons[i] in hexagons_used:
            count = count + 1
    return count
#create list of each hexagon with the amount of surrounding hexagons
unique_hexagons = pd.Series(data['start_hexagon'].unique())
list_surrounding_hexagons = unique_hexagons.apply(lambda x: Count_Surrounding_Hexagons(x,unique_hexagons))
hexagon_data = pd.DataFrame(
    {'unique_hexagons':
unique_hexagons,'surrounding_hexagons':list_surrounding_hexagons
    }, columns=['unique_hexagons','surrounding_hexagons'])



#function to calculate distance between two points
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
#function to calculate the minimum distance of a point to several points
def minimum_distance_to_multiple_points(latitudes,longitudes,latitude,longitude):
    minimum_distance = 100000000
    for i in range(len(latitudes)):

        distance_i = haversine([latitudes[i],longitudes[i]],[latitude,longitude])
        if  distance_i < minimum_distance:
            minimum_distance = distance_i
            minimum_distance_coordinates = [latitudes[i],longitudes[i]]
    return minimum_distance, minimum_distance_coordinates

Center_Amsterdam = [52.372319, 4.889525]
#functions to create columns of distances to certain landmark
hexagon_data['dist_from_center'] = hexagon_data.apply(lambda row: haversine(h3.h3_to_geo(row['unique_hexagons']),Center_Amsterdam),axis=1)


hexagon_data['dist_closest_trainstation'] = hexagon_data.apply(lambda row: minimum_distance_to_multiple_points(
                                        trainstations_dwh_ams['latitude'],trainstations_dwh_ams['longitude']
                                        ,h3.h3_to_geo(row['unique_hexagons'])[0],h3.h3_to_geo(row['unique_hexagons'])[1])[0], axis=1)
hexagon_data['distance_closest_uni_hbo'] = hexagon_data.apply(lambda row: minimum_distance_to_multiple_points(
                                       uni_hbo_dwh_ams['latitude'],uni_hbo_dwh_ams['longitude']
                                        ,h3.h3_to_geo(row['unique_hexagons'])[0],h3.h3_to_geo(row['unique_hexagons'])[1])[0], axis=1)
hexagon_data['distance_closest_metro'] = hexagon_data.apply(lambda row: minimum_distance_to_multiple_points(
                                       metro_dwh_ams['latitude'],metro_dwh_ams['longitude']
                                        ,h3.h3_to_geo(row['unique_hexagons'])[0],h3.h3_to_geo(row['unique_hexagons'])[1])[0], axis=1)
hexagon_data =hexagon_data.rename(columns={'unique_hexagons': 'start_hexagon'})
#merge hexagon data with data
data = data.merge(hexagon_data, on='start_hexagon', how='left')




#1_ring hexagon data
#function to generate data of the 6 surrounding hexagons (if they exist)
def Data_Surrounding_Hexagons(hexagon,list_of_hexagons):
    hexagons_used = list(list_of_hexagons)
    #list of surrounding hexagons
    surrounding_hexagons =list(h3.k_ring(hexagon,1))[1:]
    count = 0
    #initialize
    surrounding_data_list_1  = pd.DataFrame()
    surrounding_data_list_2  = pd.DataFrame()
    surrounding_data_list_3  = pd.DataFrame()
    a = pd.Series()
    b = pd.Series()
    c = pd.Series()
    #loop over the hexagons
    for i in range(0,6):
        if surrounding_hexagons[i] in hexagons_used:
            #if count == 0 initialize dataframe otherwise append dataframe
            if count == 0:
                surrounding_data_list = data[data.start_hexagon == surrounding_hexagons[i]].drop(data[data.start_hexagon == surrounding_hexagons[i]].index[1:])
                surrounding_data_list_1 = surrounding_data_list[data.columns[16:25]]
                surrounding_data_list_2 = surrounding_data_list[data.columns[25:61]]
                surrounding_data_list_3 = surrounding_data_list[data.columns[61:]]
            else:             
                surrounding_data_list_1 = pd.concat([surrounding_data_list_1, 
                                                     data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].drop(data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].index[1:])[data.columns[16:25]]])
                surrounding_data_list_2 = pd.concat([surrounding_data_list_2, 
                                                     data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].drop(data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].index[1:])[data.columns[25:61]]])
                surrounding_data_list_3 = pd.concat([surrounding_data_list_3, 
                                                     data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].drop(data[data.start_hexagon == 
                                                    surrounding_hexagons[i]].index[1:])[data.columns[61:]]])
            count = count + 1
    #here we take the max, min or mean of the data of the surrounding hexagons depending on the surrounding hexagons
    if not surrounding_data_list_1.empty:
        a = surrounding_data_list_1.max()
    if not surrounding_data_list_2.empty:
        b = surrounding_data_list_2.mean()
    if not surrounding_data_list_3.empty:
        c = surrounding_data_list_3.min()
    #concatenate lists
    a = pd.Series([hexagon], index=['start_hexagon']).append(a.append(b.append(c)))
    return pd.DataFrame(a).transpose()

#create the neighbouring hexagon data for all hexagons in the data set
data_neighbouring_hexagons = pd.DataFrame()
unique_hexagons = pd.Series(data['start_hexagon'].unique())
for i in range(len(unique_hexagons)):
    data_neighbouring_hexagons = data_neighbouring_hexagons.append(Data_Surrounding_Hexagons(unique_hexagons[i],unique_hexagons))
#give a new column name that differentiates the columns of the neighbouring hexagon data of the hexagon data
data_neighbouring_hexagons.columns = [col + '_neigh_hex' for col in data_neighbouring_hexagons.columns]
data_neighbouring_hexagons =data_neighbouring_hexagons.rename(columns={'start_hexagon_neigh_hex': 'start_hexagon'})
#merge the data
data = data.merge(data_neighbouring_hexagons, on='start_hexagon', how='left')

# %%add fleetsize data (STILL to do)
#first sum over datapoints at the same timestamp
fleetsize =operationality.groupby(['calculated_at']).agg({'count' : ['sum']}).reset_index()
fleetsize.columns = [col[0] for col in fleetsize.columns.values]
#then group by day as it differs throughout the day
fleetsize_day = fleetsize.groupby([pd.Grouper(key='calculated_at',freq='d')]).agg({'count' : ['mean']})
fleetsize_day.columns = [col[0] for col in fleetsize_day.columns.values]


covid = covid.drop([covid.columns[0],covid.columns[3],covid.columns[4]],1)
covid['Date'] = covid.apply(lambda row: datetime.strptime(row['Date'], '%d/%m/%Y %H:%M:%S'),axis=1)
covid['Date'][0].strftime('%d/%m/%Y %H:%M:%S')
aaa = data.merge(covid, on='start_hexagon', how='left')


# %% Visualizations NO NEED TO BE CHECKED
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

#loading in environmental data
#52.290819,4.786194,52.410362,5.038948
#http://bboxfinder.com/#52.290819,4.786194,52.410362,5.038948
# EPSG: 28992 Amersfoort / RD New
bbox = (113995.7727,478262.1749,131303.8192,491450.0742)

bbox = (532796.6789,6852881.0401,560933.1255,6874666.9972)
bbox = (78926.7613,429312.5322,104500.6136,445566.2458)
bbox = (334916.4533,6735309.9095,804545.5551,6956060.0472)

buurten = gpd.read_file(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\buurt_2020_v1.shp")
buurten.crs = "EPSG:4326"
buurten_rotterdam= buurten.to_crs(epsg=4326)
buurten_rotterdam_web = buurten_rotterdam.to_crs(epsg=3857)

buurten_rotterdam_web.save()
# buurten_rotterdam.to_pickle('BUURTEN_CBS2017_RDM_')

# CBS Squares
data500 = gpd.read_file("CBS_vk500_2020_v1.shp",bbox=bbox)
squares500_rotterdam = data500.to_crs(epsg=4326)
squares500_rotterdam_web = squares500_rotterdam.to_crs(epsg=3857)
# squares500_rotterdam.to_pickle('SQUARES_CBS2016_RDM_')

wijkbuurtkaart = gpd.read_file('Documents\GitHub\hello-world2\WijkBuurtkaart_2020_v1.gpkg')
wijkbuurtkaart_ams = wijkbuurtkaart[wijkbuurtkaart['gemeentenaam'] == 'Amsterdam']
wijkbuurtkaart_ams.crs = "EPSG:4326"
wijkbuurtkaart_ams= wijkbuurtkaart_ams.to_crs(epsg=4326)
wijkbuurtkaart_ams_web = wijkbuurtkaart_ams.to_crs(epsg=3857)

wijkbuurtkaart10 = wijkbuurtkaart.head(10)
wijkbuurtkaart_ams = wijkbuurtkaart[wijkbuurtkaart['gemeentenaam'] == 'Amsterdam']

cbs_vk500_2020 = gpd.read_file('Documents\GitHub\hello-world2\cbs_vk500_2020.gpkg')
cbs_vk500_2020_head = cbs_vk500_2020.head(1)
cbs_vk500_2020_1000 = cbs_vk500_2020.head(1000)
cbs_vk500_2020.crs = "EPSG:4326"
cbs_vk500_2020 = cbs_vk500_2020.to_crs(epsg=4326)

cbs_vk500_2020_head['geometry'] = cbs_vk500_2020_head['geometry'].apply(json.loads)

#get current working directory
os.getcwd()


data3 = data2.to_crs(epsg=4326)

a = cbs_vk500_2020_head['geometry'].convex_hull
b = a.__geo_interface__
c = b['features'][1]
d = c['geometry']
cbs_vk500_2020_head['geometry'].to_json()
list(h3.polyfill(a, 9))
b = a['features'][1]
c = b['geometry']
d = cbs_vk500_2020_head['geometry'].convex_hull
e=d.__geo_interface__
d.to_json()
def plot_zips(dataset,variable,title):

    vmin, vmax = 0, dataset[variable].max() 
    
    # create figure and axes for Matplotlib
    fig, ax = plt.subplots(1, figsize=(30, 10))

    dataset.plot(column=variable,cmap='PiYG', linewidth=0.8, ax=ax, edgecolor='0.8', alpha = 0.6) 
    ctx.add_basemap(ax, source=ctx.providers.CartoDB.Voyager)
    
    ax.axis('off')
    ax.set_title(title, fontdict={'fontsize': '20', 'fontweight' : '3'})

    # Create colorbar as a legend
    sm = plt.cm.ScalarMappable(cmap='PiYG', norm=plt.Normalize(vmin=vmin, vmax=vmax))

    # empty array for the data range
    sm._A = []

    # add the colorbar to the figure
    cbar = fig.colorbar(sm)
    plt.show()
    plt.savefig('C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/testfigure.png')
    
input_graph = CBS_vk500_2020_v1

plot_zips(input_graph,'INWONER','percentage_gehuwd households per neighborhood')

sf = gpd.read_file('C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/woonbrt10_region.shp')
with shapefile.Reader(r'C:/Users/Martijn Oerlemans/Documents/GitHub/hello-world2/woonbrt10_region.shp') as shp:
    print(shp)
    
s = sf.shape
geoj = s.__geo_interface__

geoj["type"]
polygons_used['st_asgeojson'] = polygons_used['st_asgeojson'].apply(json.loads)
geom = [shape(i) for i in polygons_used['st_asgeojson']]
polygons_used = gpd.GeoDataFrame(polygons_used,geometry=geom)
plot_zips(input_graph,'percentage_gehuwd','percentage_gehuwd households per neighborhood')

input_graph = cbs_vk500_2020

# %%
# save the resulting dataframe to gcp
filename = 'forecasting_demand_df_{}'.format(datetime.now().strftime("%d:%m:%Y-%H:%M:%S"))
folder = 'Real-time forecasting demand/'
write_pickle_to_bucket(data,'felyx-ai-machine-learning',folder,filename)
