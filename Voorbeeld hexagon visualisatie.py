# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:07:30 2021

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
#set up database connection
os.environ["GOOGLE_APPLICATION_CREDENTIALS"]="C:/Users/Martijn Oerlemans/Desktop/felyx-ai-machine-learning-88e8adf81f5b.json"
cnx = create_engine(access_secret_version('felyx-ai-machine-learning','DATABASE_URL_DATAWAREHOUSE', "latest"))
#load in data
closest_vehicle_on_app_resume = pd.read_sql_query('''SELECT received_at, reservation_created, reservation_tripstarted, user_id,  user_longitude, user_latitude
                                       FROM ios_production.closest_vehicle_on_app_resume
                                       WHERE user_latitude IS NOT NULL
                                        AND user_latitude IS NOT NULL
                                       AND DATE_TRUNC('day',received_at) > '2021-02-28'
                                       AND DATE_TRUNC('day',received_at) < '2021-08-02'
                                       ''',cnx)
#set size of hexagon                               
service_area_resolution = 9
#add the hexagon column based on latitudes and longitudes
closest_vehicle_on_app_resume['hexagon'] = closest_vehicle_on_app_resume.apply(lambda row: 
                                    h3.geo_to_h3(lat=row['user_latitude'],
                                    lng=row['user_longitude'], 
                                    resolution=service_area_resolution), axis=1)
#function to plot the hexagons, with as input the hexagons, their color and whether you want to 
#create a new map or add the hexagons to an already existing map  
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
        #here you can determine what the opacity is and the thickness of the borders of the hexagon
        my_PolyLine=folium.PolyLine(locations=polyline,weight=0.5,color=color,fill=True,fill_color = color,fill_opacity = 0.6)
        m.add_child(my_PolyLine)
    return m

#Determining the color range
app_opening_count =closest_vehicle_on_app_resume['hexagon'].value_counts()
maximum_opening_hexagon = max(app_opening_count)
#creating the color range
def colorFader(c1,c2,mix): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)
#colors wanted
c1 = 'green'
c2='red'
#create the hexagon map and loop over the different app_opening_counts so that
#each hexagon has a color corresponding to their app_opening_counts value
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
#create the legenda
colormap = branca.colormap.LinearColormap(colors=['green','red'],vmin=0,vmax=maximum_opening_hexagon)

colormap = colormap.to_step(index=[0, maximum_opening_hexagon* (1/10), maximum_opening_hexagon* (2/10)
                                   , maximum_opening_hexagon* (3/10), maximum_opening_hexagon* (4/10)
                                   , maximum_opening_hexagon* (5/10), maximum_opening_hexagon* (6/10)
                                   , maximum_opening_hexagon* (7/10), maximum_opening_hexagon* (8/10),
                                   maximum_opening_hexagon* (9/10), maximum_opening_hexagon])
colormap.caption = 'App openings last week per hexagon (service_area)'
colormap.add_to(m)
#Add a marker where we show the total app openings
folium.Marker(location=[51.920731882528166, 4.470521006248073],popup='Total app openings :' 
              + str(len(closest_vehicle_on_app_resume)),).add_to(m)
#Save the map
m.save(r"C:\Users\Martijn Oerlemans\Documents\GitHub\hello-world2\app_openings_last_week1.html")