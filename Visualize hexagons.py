# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 16:17:56 2021

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
#load data from sql database
from felyx_gcp_utils.get_gcp_secrets import access_secret_version
import json
from shapely.geometry import shape
#import geoplot
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry import Polygon
from shapely.geometry import MultiPoint
import geopandas as gpd
os.environ["GOOGLE_APPLICATION_CREDENTIALS"]="C:/Users/Martijn Oerlemans/Desktop/felyx-ai-machine-learning-88e8adf81f5b.json"
cnx = create_engine(access_secret_version('felyx-ai-machine-learning','DATABASE_URL_DATAWAREHOUSE', "latest"))

closest_vehicle_on_app_resume_monday = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-03-02'
                                       ''',cnx)
closest_vehicle_on_app_resume = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       ''',cnx)
closest_vehicle_on_app_resume_morning = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) < 10
                                       ''',cnx)
closest_vehicle_on_app_resume_noon = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) > 9
                                       AND DATE_PART('hour',received_at) < 15                                       
                                       ''',cnx)
closest_vehicle_on_app_resume_night = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) > 14                                      
                                       ''',cnx)
reservation_park_mode_enabled = pd.read_sql_query('''SELECT received_at, user_id,  user_longitude, user_latitude
                                       FROM ios_production.reservation_park_mode_enabled 
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       ''',cnx)
reservation_park_mode_enabled_morning = pd.read_sql_query('''SELECT received_at, user_id,  user_longitude, user_latitude
                                       FROM ios_production.reservation_park_mode_enabled
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) < 10
                                       ''',cnx)
reservation_park_mode_enabled_noon = pd.read_sql_query('''SELECT received_at, user_id,  user_longitude, user_latitude
                                       FROM ios_production.reservation_park_mode_enabled 
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) > 9
                                       AND DATE_PART('hour',received_at) < 15                                       
                                       ''',cnx)
reservation_park_mode_enabled_night = pd.read_sql_query('''SELECT received_at, user_id,  user_longitude, user_latitude
                                       FROM ios_production.reservation_park_mode_enabled 
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       AND DATE_PART('hour',received_at) > 14                                      
                                       ''',cnx)
                                       
service_area = pd.read_sql_query('''SELECT *
                                FROM public.service_area
                                WHERE location_id =1                                     
                                       ''',cnx)
                                       
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

polygons = pd.read_sql_query('''SELECT ST_AsGeoJSON(ST_Transform(wgs84_polygon,4326)), id, default_name 
                             from public.service_area
                             WHERE location_id=1''',cnx)
closest_vehicle_on_app_resume_monday['hexagon'] = closest_vehicle_on_app_resume_monday.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
closest_vehicle_on_app_resume['hexagon'] = closest_vehicle_on_app_resume.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
closest_vehicle_on_app_resume_morning['hexagon'] = closest_vehicle_on_app_resume_morning.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)  
closest_vehicle_on_app_resume_noon['hexagon'] = closest_vehicle_on_app_resume_noon.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)  
closest_vehicle_on_app_resume_night['hexagon'] = closest_vehicle_on_app_resume_night.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
reservation_park_mode_enabled['hexagon'] = reservation_park_mode_enabled.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
reservation_park_mode_enabled_morning['hexagon'] = reservation_park_mode_enabled_morning.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)   
reservation_park_mode_enabled_noon['hexagon'] = reservation_park_mode_enabled_noon.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)   
reservation_park_mode_enabled_night['hexagon'] = reservation_park_mode_enabled_night.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
rides_unique_polygons = rides.service_area_start.unique()

polygons_used = polygons[polygons['default_name'].isin(rides_unique_polygons)]
polygons_used_unique = polygons_used.default_name.unique()
polygons_used['st_asgeojson'] = polygons_used['st_asgeojson'].apply(json.loads)
geom = [shape(i) for i in polygons_used['st_asgeojson']]
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
        my_PolyLine=folium.PolyLine(locations=polyline,weight=0.5,color=color,fill=True,fill_color = color,fill_opacity = 0.6)
        m.add_child(my_PolyLine)
    return m

app_opening_count_monday =closest_vehicle_on_app_resume_monday['hexagon'].value_counts()
maximum_opening_hexagon_monday = max(app_opening_count_monday)
app_opening_count =closest_vehicle_on_app_resume['hexagon'].value_counts()
maximum_opening_hexagon = max(app_opening_count)
app_opening_count_morning =closest_vehicle_on_app_resume_morning['hexagon'].value_counts()
maximum_opening_hexagon_morning = max(app_opening_count_morning)
app_opening_count_noon =closest_vehicle_on_app_resume_noon['hexagon'].value_counts()
maximum_opening_hexagon_noon = max(app_opening_count_noon)
app_opening_count_night =closest_vehicle_on_app_resume_night['hexagon'].value_counts()
maximum_opening_hexagon_night = max(app_opening_count_night)

reservation_park_mode_enabled_count =reservation_park_mode_enabled['hexagon'].value_counts()
maximum_reservation_park_mode_enabled_hexagon = max(reservation_park_mode_enabled_count)

def colorFader(c1,c2,mix): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

#whole wheek app openings
c1 = 'green'
c2='red'

for i in range(len(app_opening_count.unique())):
    if i == 0:
        m = visualize_hexagons(app_opening_count.index[app_opening_count 
                               == app_opening_count.unique()[0]],
                               colorFader(c1,c2,
                                app_opening_count.unique()[0]/maximum_opening_hexagon) 
                               ,None)
    else:
        m = visualize_hexagons(app_opening_count.index[app_opening_count 
                               == app_opening_count.unique()[i]],
                               colorFader(c1,c2,
                                app_opening_count.unique()[i]/maximum_opening_hexagon) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_opening_hexagon)

colormap = colormap.to_step(index=[0, maximum_opening_hexagon* (1/10), maximum_opening_hexagon* (2/10)
                                   , maximum_opening_hexagon* (3/10), maximum_opening_hexagon* (4/10)
                                   , maximum_opening_hexagon* (5/10), maximum_opening_hexagon* (6/10)
                                   , maximum_opening_hexagon* (7/10), maximum_opening_hexagon* (8/10),
                                   maximum_opening_hexagon* (9/10), maximum_opening_hexagon])
colormap.caption = 'App openings last week per hexagon (service_area)'
colormap.add_to(m)
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total app openings :' 
              + str(len(closest_vehicle_on_app_resume)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\app_openings_last_week.html")
#morning app openings
c1 = 'green'
c2='red'

for i in range(len(app_opening_count_morning.unique())):
    if i == 0:
        m = visualize_hexagons(app_opening_count_morning.index[app_opening_count_morning 
                               == app_opening_count_morning.unique()[0]],
                               colorFader(c1,c2,
                                app_opening_count_morning.unique()[0]/maximum_opening_hexagon_morning) 
                               ,None)
    else:
        m = visualize_hexagons(app_opening_count_morning.index[app_opening_count_morning 
                               == app_opening_count_morning.unique()[i]],
                               colorFader(c1,c2,
                                app_opening_count_morning.unique()[i]/maximum_opening_hexagon_morning) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_opening_hexagon_morning)

colormap = colormap.to_step(index=[0, maximum_opening_hexagon_morning* (1/10), maximum_opening_hexagon_morning* (2/10)
                                   , maximum_opening_hexagon_morning* (3/10), maximum_opening_hexagon_morning* (4/10)
                                   , maximum_opening_hexagon_morning* (5/10), maximum_opening_hexagon_morning* (6/10)
                                   , maximum_opening_hexagon_morning* (7/10), maximum_opening_hexagon_morning* (8/10),
                                   maximum_opening_hexagon_morning* (9/10), maximum_opening_hexagon_morning])
colormap.caption = 'App openings last week in the morning (until 11 am) per hexagon (service_area)'
colormap.add_to(m)
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total app openings in the morning:' 
              + str(len(closest_vehicle_on_app_resume_morning)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\app_openings_last_week_morning.html")
#noon app openings
c1 = 'green'
c2='red'

for i in range(len(app_opening_count_noon.unique())):
    if i == 0:
        m = visualize_hexagons(app_opening_count_noon.index[app_opening_count_noon 
                               == app_opening_count_noon.unique()[0]],
                               colorFader(c1,c2,
                                app_opening_count_noon.unique()[0]/maximum_opening_hexagon_noon) 
                               ,None)
    else:
        m = visualize_hexagons(app_opening_count_noon.index[app_opening_count_noon 
                               == app_opening_count_noon.unique()[i]],
                               colorFader(c1,c2,
                                app_opening_count_noon.unique()[i]/maximum_opening_hexagon_noon) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_opening_hexagon_noon)

colormap = colormap.to_step(index=[0, maximum_opening_hexagon_noon* (1/10), maximum_opening_hexagon_noon* (2/10)
                                   , maximum_opening_hexagon_noon* (3/10), maximum_opening_hexagon_noon* (4/10)
                                   , maximum_opening_hexagon_noon* (5/10), maximum_opening_hexagon_noon* (6/10)
                                   , maximum_opening_hexagon_noon* (7/10), maximum_opening_hexagon_noon* (8/10),
                                   maximum_opening_hexagon_noon* (9/10), maximum_opening_hexagon_noon])
colormap.caption = 'App openings last week in the noon (11 am- 3 pm) per hexagon (service_area)'
colormap.add_to(m)
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total app openings in the noon:' 
              + str(len(closest_vehicle_on_app_resume_noon)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\app_openings_last_week_noon.html")
#night app openings
c1 = 'green'
c2='red'

for i in range(len(app_opening_count_night.unique())):
    if i == 0:
        m = visualize_hexagons(app_opening_count_night.index[app_opening_count_night
                               == app_opening_count_night.unique()[0]],
                               colorFader(c1,c2,
                                app_opening_count_night.unique()[0]/maximum_opening_hexagon_night) 
                               ,None)
    else:
        m = visualize_hexagons(app_opening_count_night.index[app_opening_count_night 
                               == app_opening_count_night.unique()[i]],
                               colorFader(c1,c2,
                                app_opening_count_night.unique()[i]/maximum_opening_hexagon_night) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_opening_hexagon_night)

colormap = colormap.to_step(index=[0, maximum_opening_hexagon_night* (1/10), maximum_opening_hexagon_night* (2/10)
                                   , maximum_opening_hexagon_night* (3/10), maximum_opening_hexagon_night* (4/10)
                                   , maximum_opening_hexagon_night* (5/10), maximum_opening_hexagon_night* (6/10)
                                   , maximum_opening_hexagon_night* (7/10), maximum_opening_hexagon_night* (8/10),
                                   maximum_opening_hexagon_night* (9/10), maximum_opening_hexagon_night])
colormap.caption = 'App openings last week in the night (> 3 pm) per hexagon (service_area)'
colormap.add_to(m)
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total app openings in the night:' 
              + str(len(closest_vehicle_on_app_resume_night)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\app_openings_last_week_night.html")

#whole wheek reservation park mode enabled
c1 = 'green'
c2='red'

for i in range(len(reservation_park_mode_enabled_count.unique())):
    if i == 0:
        m = visualize_hexagons(reservation_park_mode_enabled_count.index[reservation_park_mode_enabled_count
                               == reservation_park_mode_enabled_count.unique()[0]],
                               colorFader(c1,c2,
                                reservation_park_mode_enabled_count.unique()[0]/maximum_reservation_park_mode_enabled_hexagon) 
                               ,None)
    else:
        m = visualize_hexagons(reservation_park_mode_enabled_count.index[reservation_park_mode_enabled_count 
                               == reservation_park_mode_enabled_count.unique()[i]],
                               colorFader(c1,c2,
                                reservation_park_mode_enabled_count.unique()[i]/maximum_reservation_park_mode_enabled_hexagon) 
                               ,m)

colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_reservation_park_mode_enabled_hexagon)

colormap = colormap.to_step(index=[0, maximum_reservation_park_mode_enabled_hexagon* (1/10), maximum_reservation_park_mode_enabled_hexagon* (2/10)
                                   , maximum_reservation_park_mode_enabled_hexagon* (3/10), maximum_reservation_park_mode_enabled_hexagon* (4/10)
                                   , maximum_reservation_park_mode_enabled_hexagon* (5/10), maximum_reservation_park_mode_enabled_hexagon* (6/10)
                                   , maximum_reservation_park_mode_enabled_hexagon* (7/10), maximum_reservation_park_mode_enabled_hexagon* (8/10),
                                   maximum_reservation_park_mode_enabled_hexagon* (9/10), maximum_reservation_park_mode_enabled_hexagon])
colormap.caption = 'Park mode enabled last week per hexagon (service_area)'
colormap.add_to(m)
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total reservations park mode:' 
              + str(len(reservation_park_mode_enabled)),).add_to(m)
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\reservation_enabled_last_week.html")
