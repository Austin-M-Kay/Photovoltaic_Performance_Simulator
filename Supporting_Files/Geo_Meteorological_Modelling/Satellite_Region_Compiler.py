# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 18:18:31 2026

@author: austi

#==============================================================================
# Compile All Information About Each Satellite Region Into a JSON Object
#==============================================================================

Using the Satellite Region names as keys, store the following attributes:
    
    - Available years
    - Model name
    - Source
    - Temporal resolution options
    - Spatial resolution
    - URL start
    - Longitude/Latitudes forming the perimeter
    
"""

#------------------------------------------------------------------------------
# Preliminary Imports
#------------------------------------------------------------------------------

# Will need to install this module if not present:
# pip install shapely

from json import dump
from numpy import arange, array
from pandas import read_excel
import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Load Satellite Overview Data
#------------------------------------------------------------------------------

Satellite_Overview_DataFrame = read_excel( "Satellite_Overview.xlsx" )

Satellite_Names = list( Satellite_Overview_DataFrame[ "Region" ] )

Output = {}

plt.figure( dpi = 200 )

for j, Satellite in enumerate( Satellite_Names ):

    Row = Satellite_Overview_DataFrame.iloc[ j ]

    #--------------------------------------------------------------------------
    # Store Row Parameters
    #--------------------------------------------------------------------------    
    
    Output_j = { Key : str( Row[ Key ] ) for Key in list( Row.keys() ) 
                
                if Key not in [ "Years", "Region", "Temporal Resolution [min]" ] }
    
    #--------------------------------------------------------------------------
    # Store Temporal Resolutions
    #--------------------------------------------------------------------------    
    
    Output_j[ "Temporal Resolution [min]" ] = str( Row[ "Temporal Resolution [min]" ] ).split( ", " )
    
    #--------------------------------------------------------------------------
    # Store Years
    #--------------------------------------------------------------------------    
    
    Years = Row[ "Years" ].split( "-" )
        
    Years = arange( 
        int( min( Years ) ), 
        int( max( Years ) ) + 1 )
    
    Years = [ str( Year ) for Year in Years ]
    
    Output_j[ "Years" ] = Years    
        
    Output[ Satellite ] = Output_j 

    #--------------------------------------------------------------------------
    # Define Coordinates
    #--------------------------------------------------------------------------    
        
    Perimeter_Data = read_excel( Satellite + '.xlsx' ) 
    
    Longitudes = list( Perimeter_Data[ "Longitude" ] )
    Latitudes = list( Perimeter_Data[ "Latitude" ] )
        
    # Create tuples of coordinates
    
    Coordinates = list( array( [ 
        Longitudes, 
        Latitudes ] ).transpose().tolist() )
    
    #--------------------------------------------------------------------------
    # Store Coordinates and Plot
    #--------------------------------------------------------------------------    
    
    Output_j[ "Coordinates" ] = Coordinates    

    plt.plot( 
        Longitudes + [ Longitudes[ 0 ] ],
        Latitudes + [ Latitudes[ 0 ] ], 
        # label = Satellite 
        )    
    
    #--------------------------------------------------------------------------
    # Some Satellites Cross the -180/+180 Line, Determine, Store and Plot These 
    #--------------------------------------------------------------------------    
    
    if "Latitude Continued" in Perimeter_Data.columns:
        
        # Identify
        
        Data_Continued = Perimeter_Data[ [ "Longitude Continued", "Latitude Continued" ] ].dropna()
        
        Longitudes_Continued = list( Data_Continued[ "Longitude Continued" ] )
        Latitudes_Continued  = list( Data_Continued[ "Latitude Continued" ] )
        
        # Define Coordinates
        
        Coordinates_Continued = list( array( [ 
            Longitudes_Continued, 
            Latitudes_Continued ] ).transpose().tolist() )
        
        # Store

        Output_j[ "Coordinates Continued" ] = Coordinates_Continued 

    #--------------------------------------------------------------------------
    # Store Output
    #--------------------------------------------------------------------------
    
    Output[ Satellite ] = Output_j 


plt.xlabel( "Longitude" )
plt.ylabel( "Latitude" )
plt.show()  

#------------------------------------------------------------------------------
# Store Output
#------------------------------------------------------------------------------
        
with open( "Satellite_Information.json", "w" ) as File:
    dump( Output, File, indent = 2 )
    File.close()