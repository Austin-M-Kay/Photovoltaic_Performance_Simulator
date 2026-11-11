# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 07:28:17 2026

@author: austi
"""

# pip install shapely

from pandas import read_excel, notnull
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from numpy import array
import matplotlib.pyplot as plt

Satellite_Names = [
    "Americas",
    "GOES East - Extended Americas",
    "GOES West - Extended Americas",
    "Arctic",
    "Extended CONUS - USA and Mexico",
    "METEOSAT IODC - Europe, Africa, and Asia",
    "METEOSAT Prime Meridian - Africa and Europe",
    "Himawari 2km Resolution - East-Asia, Australia, and Pacific",
    "Himawari 4km Resolution - East-Asia, Australia, and Pacific",
    "SUNY - South Asia" ]

Satellite_Polygons = {}

plt.figure( dpi = 200 )
for Satellite in Satellite_Names:
    
    Data = read_excel( Satellite + '.xlsx' ) 
    
    Longitudes = list( Data[ "Longitude" ] )
    Latitudes = list( Data[ "Latitude" ] )
        
    # If continuous over -180 to 180 longitude range
    
    Satellite_Polygons[ Satellite ] = Polygon( 
        array( [ 
            Longitudes, 
            Latitudes ] ).transpose().tolist() )
        
    plt.plot( 
        Longitudes + [ Longitudes[ 0 ] ],
        Latitudes + [ Latitudes[ 0 ] ], 
        # label = Satellite 
        )    
    
    if "Latitude Continued" in Data.columns:
        Data_Continued = Data[ [ "Longitude Continued", "Latitude Continued" ] ].dropna()
        
        Longitudes_Continued = list( Data_Continued[ "Longitude Continued" ] )
        Latitudes_Continued  = list( Data_Continued[ "Latitude Continued" ] )

        Satellite_Polygons[ Satellite + '_Continued' ] = Polygon( 
            array( [ 
                Longitudes_Continued, 
                Latitudes_Continued ] ).transpose().tolist() )
        
        plt.plot( 
            Longitudes_Continued + [ Longitudes_Continued[ 0 ] ],
            Latitudes_Continued + [ Latitudes_Continued[ 0 ] ], 
            # label = Satellite 
            )    

# plt.legend()
plt.xlabel( "Longitude" )
plt.ylabel( "Latitude" )
plt.show()

def Satellite_Assigner(
    Latitude,
    Longitude,
    Satellite_Polygons = Satellite_Polygons,
    ):
    
    Location_Point = Point( Longitude, Latitude )
    
    Assigned_Satellites = []
    
    for Satellite in list( Satellite_Polygons.keys() ):
        print( Satellite, Satellite_Polygons[ Satellite ].contains( Location_Point ) )

        if Satellite_Polygons[ Satellite ].contains( Location_Point ):
            Assigned_Satellites.append( Satellite.replace( "_Continued" , "" ) )
    
    return Assigned_Satellites


Satellite_Assigner( 20, -155)


#%%

import requests
from pandas import read_csv 

# METEOSAT IODC
# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/msg-iodc-download.csv?"

# METEOSAT PRime Mereidian

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/nsrdb-msg-v1-0-0-download.csv?"
# Lat = 2.2 
# Lon = 0 
# Year = 2018 

# Polar 

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/nsrdb-polar-v4-0-0-download.csv?"
# Lat = 80.4 
# Lon = 80.5 
# Year = 2020 

# Asia - 2011-2015

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/himawari7-download.csv?"
# Lat = 36 
# Lon = 138 
# Year = 2011 

# Asia - 2016-2020

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/himawari-download.csv?"
# Lat = 36 
# Lon = 138 
# Year = 2020 

# US CONUS - 2018-2025

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/nsrdb-GOES-conus-v4-0-0-download.csv?"
# Lat = 45 
# Lon = -70
# Year = 2025

# SUNY - India 2000-2014

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/suny-india-download.csv?"
# Lat = 30
# Lon = 75
# Year = 2000

# US CONUS - 1998-2025

# url = "https://developer.nlr.gov/api/nsrdb/v2/solar/nsrdb-GOES-aggregated-v4-0-0-download.csv?"
# Lat = 0
# Lon = -70
# Year = 1998

# US Full Disc - 2018-2025

url = "https://developer.nlr.gov/api/nsrdb/v2/solar/nsrdb-GOES-full-disc-v4-0-0-download.csv?"
Lat = 57
Lon = -156
Year = 2018


#--------------------------------

API_Key = "w7n6NcPPqTzPh7qJFFcj99JbMz8rmpWJ6LccNFjd"

url += "api_key=" + API_Key

payload = [
    "names=" + str( Year ),
    "leap_day=false",
    "interval=60",
    "utc=false",
    # "full_name=Honored%2BUser"
    "email=a.m.kay@swansea.ac.uk",
    # "affiliation=NREL",
    # "mailing_list=true",
    "reason=Academic",
    # "attributes=dhi%2Cdni%2Cwind_speed%2Cair_temperature",
    "wkt=POINT(" + str( Lon ) + "%20" + str( Lat ) + ")",
    # "wkt=MULTIPOINT(-106.22%2032.9741%2C-106.18%2032.9741%2C-106.1%2032.9741)"
    ]

payload = '&'.join( payload )

headers = {
    'content-type': "application/x-www-form-urlencoded",
    'cache-control': "no-cache"
}

response = requests.request("POST", url, data=payload, headers=headers)

Data = [ Row.split( ',' ) for Row in response.text.split( '\n' ) ]

Meta_Data = { Head : Data[ 1 ][ j ] for j, Head in enumerate( Data[ 0 ] ) }

Columns = Data[ 2 ]
Data = Data[ 3: ]


from pandas import DataFrame

DF = DataFrame( Data, columns = Columns )

DF.to_clipboard( index = False )

print( DF )























