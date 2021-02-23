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