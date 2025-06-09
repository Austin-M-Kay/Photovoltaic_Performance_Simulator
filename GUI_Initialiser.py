# -*- coding: utf-8 -*-
"""
Created on Fri May 23 22:23:57 2025

@author: Austin M. Kay
"""

from os import getcwd, path

Filepath = path.join( getcwd(), "Initialise_GUI.txt" )

with open( Filepath, "w" ) as File:
    
    File.close()
