# -*- coding: utf-8 -*-
"""
Created on Wed Apr 30 23:15:31 2025

@author: austi
"""

from os import getcwd, path, mkdir
from pandas import read_excel

CWD = getcwd()

#Filename for library here:
    
DF = read_excel( "ADD_Optical_Constants.xlsx" )

Column_Headers = list( DF.columns )

Material_Types = [ Material.replace( '_lambda', '' ) for Material in Column_Headers ]
Material_Types = [ Material.replace( '_n', '' ) for Material in Material_Types ]
Material_Types = [ Material.replace( '_k', '' ) for Material in Material_Types ]

Material_Set = list( set( Material_Types ) )

Counts = { Material : 0 for Material in Material_Set }

for Material in Material_Types:
    
    Counts[ Material ] += 1
    
for Material in Material_Set:
    
    if Counts[ Material ] < 3:
        
        Material_Set.remove( Material )
        
        print( Material )
        
    if Counts[ Material ] > 3:
        
        print( Material )
        
# 'ITOsorizon_alpha': 1,        
# 'PCDTBTPCBM21': 1,

# Create Output File

Output_Path = path.join( CWD, "Separted_Constants" )

mkdir( Output_Path )

for Material in Material_Set:

    Output_n = DF[ [ Material + "_lambda", Material + '_n' ] ].dropna()
    Output_k = DF[ [ Material + "_lambda", Material + '_k' ] ].dropna()

    Filepath_n = path.join( Output_Path, Material + '_n.xlsx' )
    Filepath_k = path.join( Output_Path, Material + '_k.xlsx' )

    Output_n.to_excel( Filepath_n, index = False )
    Output_k.to_excel( Filepath_k, index = False )
    