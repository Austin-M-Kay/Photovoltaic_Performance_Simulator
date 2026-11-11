# -*- coding: utf-8 -*-
#%%
"""
Created on Wed Mar 12 10:22:26 2025

@author: Austin M. Kay

Email: a.m.kay.954708@swansea.ac.uk

Sustainable Advanced Materials (Sêr SAM) | Deunyddiau Uwch Cynaliadwy,
Centre for Intergrative Semiconductor Materials | Y Ganolfan Deunyddiau Lled-ddargludol Integreiddiol,
Department of Physics | Adran Ffiseg,
Swansea University | Prifysgol Abertawe,
Wales | Cymru,
SA1 8EN

"""

#==============================================================================
# Photovoltaic Performance Simulator (PV-Simulator) Version 3.0
#==============================================================================

"""
A tool for simulating the performance of photovoltaic devices under different
illumination conditions. 

From V3.0 onwards, the graphical user interface (GUI) of the tool is 
constructed in tkinter, as opposed to Jupyter (which was used in prior
versions).

When executed, the script will run the code in the following sections, 
including the "if" statements that produce case-specific outputs.

    0. Preliminary Imports
    1. File Identification and Parameter Loading
    2. Helper Functions
    3. Functions for Small-Area Thermodynamic Limit Calculations
    4. Functions for Simulating EQE Spectra
    5. Optical Modelling
    6. Drift-Diffusion Simulations
    7. Scaled-Up Device Simulations and Systems
    8. GUI
    9. "If" Statements
   
These collapsable sections are indicated by different headings, which include 
their numbers on the right-hand side.

This version of PV_Simulator has been entirely reworked from the ground up
to enhance scalability and ease-of-customisation. It is based on 
investigations, methods, models, and equations outlined in the following works:

(I)     Zarrabi, N., Sandberg, O.J., Zeiske, S., Li, W., Riley, D.B., Meredith, 
P., and Armin., A (2020). Charge-Generating Mid-Gap Trap States Define the
Thermodynamic Limit of Organic Photovoltaic Devices. Nat. Commun. 11, 5567. 
https://doi.org/10.1038/s41467-021-24202-9

(II)    Kaiser, C., Sandberg, O.J., Zarrabi, N., Li, W., Meredith, P., and 
Armin., A (2021). A Universal Urbach Rule for Disordered Organic Semiconductors
Nat. Commun. 12, 3988. 
https://doi.org/10.1038/s41467-021-24202-9

(III)   Kay, A.M., Sandberg, O.J., Zarrabi, N., Li, W., Zeiske, S., Kaiser, C., 
Meredith, P., and Armin, A. (2022), Quantifying the Excitonic Static Disorder 
in Organic Semiconductors. Adv. Funct. Mater. 32, 2113181. 
https://doi.org/10.1002/adfm.202113181
    
(IV)    Kay, A.M., Fitzsimons, M.E., Burwell, G., Meredith, P., Armin, A. and Sandberg, 
O.J. (2023), The Thermodynamic Limit of Indoor Photovoltaics Based on 
Energetically-Disordered Molecular Semiconductors. Sol. RRL, 7: 2300277.
https://doi.org/10.1002/solr.202300277

(V)     Seunarine, K., Haymoor, Z., Spence, M., Burwell, G., Kay, A.M., Meredith, P.,
Armin, A., and Carnie, M. (2024). Light Power Resource Availability for Energy
Harvesting Photovoltaics for Self-Powered IoT. J. Phys. Energy 6 (1), 015018
https://iopscience.iop.org/article/10.1088/2515-7655/ad1764/meta

(VI)    Kay, A.M., Riley, D.B., Meredith, P., Armin, A., and Sandberg, O.J. (2024), A
New Framework for Understanding Recombination-Limited Charge Extraction in 
Disordered Semiconductors. J. Phys. Chem. Lett.: 15 (16), 4416-4421
https://pubs.acs.org/doi/full/10.1021/acs.jpclett.4c00218

(VII)   Burwell, G., Zeiske, S., Caprioglio, P., Sandberg, O.J., Kay, A.M., Farrar, 
M.D., Kim, Y.R., Snaith, H.J., Meredith, P. and Armin, A. (2024), Wide-Gap 
Perovskites for Indoor Photovoltaics. Sol. RRL, 8: 2400180.
https://doi.org/10.1002/solr.202400180

(VIII)  Kay, A.M., Riley, D.B., Sandberg, O.J., Burwell, G., Meredith, P. and Armin, A.
(2024), On the Performance Limits of Agrivoltaics—From Thermodynamic to 
Geo-Meteorological Considerations. Sol. RRL, 8: 2400456. 
https://doi.org/10.1002/solr.202400456

"""

#%%===========================================================================0
# 0. Preliminary Imports
#==============================================================================

#%matplotlib inline

from datetime import datetime
from numpy import arange, arcsin, array, cos, empty, exp, floating, floor, inf, interp, integer, linalg, linspace, log, log10, logspace, meshgrid, nan, ndarray, pi, sin, sqrt, tanh, where, zeros
from os import getcwd, listdir, mkdir, path, remove, walk
from IPython.display import display, Image
from ipywidgets import BoundedFloatText, BoundedIntText, Button, Dropdown, HBox, Label, Layout, RadioButtons, VBox
from ipywidgets import Output as Output_Widget 
from json import load
from pandas import DataFrame, ExcelWriter, read_csv, read_excel
from scipy.constants import c, e, h, k
from scipy.optimize import brentq, minimize
from scipy.special import erf, erfc
from sys import exit
from time import sleep

import json
import matplotlib.pyplot as plt
import xlsxwriter

import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

# Depends on Python version:
try:
    from scipy.integrate import cumulative_simpson, simpson
except:
    from scipy.integrate import cumulative_trapezoid as cumulative_simpson # use trapz if simpson not available
    from scipy.integrate import simps as simpson
else:
    from scipy.integrate import cumulative_simpson, simpson

#%%===========================================================================1
# 1. File Identification and Parameter Loading
#==============================================================================

"""
In this section, supporting files are identified from the "Supporting_Files" 
directory. This includes the parameters and the avaialble EQE spectra, which 
are cross-referenced with PV-Simulator's Systems Database (System_Database.xlsx).

In addition, the light sources, optical constants, and more are identified.
"""

Current_Working_Directory = getcwd()

Supporting_Files_Path = path.join( Current_Working_Directory, "Supporting_Files")

#------------------------------------------------------------------------------
# Load Parameters
#------------------------------------------------------------------------------

# Parameters_Path = path.join( Current_Working_Directory, "Parameters.json" )

# with open( Parameters_Path, 'r' ) as File:
    
#     Parameters = json.load( File )
    
# File.close()

#------------------------------------------------------------------------------
# Load Systems Database
#------------------------------------------------------------------------------

Systems_Database_Path = path.join( Current_Working_Directory, "System_Database.xlsx" )

Systems_Database = read_excel( Systems_Database_Path ).to_dict()

System_Identifiers = list( Systems_Database[ "System_Ref" ].values() )

# Rearrange According to System Refs

Systems_Database = { Column : { 
    
    Systems_Database[ "System_Ref" ][ j ] : Systems_Database[ Column ][ j ] 
    
    for j in list( Systems_Database[ "System_Ref" ].keys() ) } for Column in list( Systems_Database.keys() ) }

# Remove Unneccessary Columns

Systems_Database.pop( 'System_Ref' )

#------------------------------------------------------------------------------
# Identify EQE Spectra
#------------------------------------------------------------------------------

EQE_Spectra_Folder_Path = path.join( Supporting_Files_Path, "EQE_Spectra" )

Available_EQEs = listdir( EQE_Spectra_Folder_Path )

#------------------------------------------------------------------------------
# Cross Reference EQE Spectra and Systems Database
#------------------------------------------------------------------------------

Systems_Database = { Column : { 
    
    Identifier : Systems_Database[ Column ][ Identifier ] 
    
    for Identifier in System_Identifiers if Identifier + '.xlsx' in Available_EQEs } 
    
    for Column in list( Systems_Database.keys() ) }

# Update System Identifiers

System_Identifiers = list( Systems_Database[ 'E_lower [eV]' ].keys() )

#------------------------------------------------------------------------------
# Determine Available Light Sources
#------------------------------------------------------------------------------

Light_Source_Path = path.join( Supporting_Files_Path, "Light_Sources" )

Light_Sources = listdir( Light_Source_Path )

Light_Sources = [ Source.replace( '.xlsx', '' ) for Source in Light_Sources if '.xlsx' in Source ] # remove any other types of files

#------------------------------------------------------------------------------
# Load Board Efficiency Data
#------------------------------------------------------------------------------

Board_Eff_Folder_Path = path.join( Supporting_Files_Path, "Energy_Harvesting_Board_Efficiencies" )

I_100uA_Data = read_excel( path.join( Board_Eff_Folder_Path, "0.1mA.xlsx" ) )
I_1mA_Data = read_excel( path.join( Board_Eff_Folder_Path, "1mA.xlsx" ) )
I_10mA_Data = read_excel( path.join( Board_Eff_Folder_Path, "10mA.xlsx" ) )

#------------------------------------------------------------------------------
# Determine Available Simulation Modes
#------------------------------------------------------------------------------

Simulation_Modes_Path = path.join( Supporting_Files_Path, "GUI_Files" )
Simulation_Modes_Path = path.join( Supporting_Files_Path, "Simulation_Modes.xlsx" ) 

Simulation_Modes = read_excel( Simulation_Modes_Path )
Simulation_Modes = list( Simulation_Modes[ "Mode Name" ] )

#------------------------------------------------------------------------------
# Determine Available Optical Constants
#------------------------------------------------------------------------------

Optical_Constants_Folder_Path = path.join( Supporting_Files_Path , 'Optical_Constants' )

Optical_File_List = []
Folder_Names = []

for Root, Dirs, Files in walk( Optical_Constants_Folder_Path ):
    Folder_Names.append( Dirs )
    for File in Files:
        Trimmed_Filename = File.replace( '.xlsx' , '' )
        Optical_File_List.append( Trimmed_Filename ) 

n_File_List = [ Entry.replace( '_n' , '' ) for Entry in Optical_File_List if '_n' in Entry ]
k_File_List = [ Entry.replace( '_k' , '' ) for Entry in Optical_File_List if '_k' in Entry ]

Materials = [ Entry for Entry in n_File_List if Entry in k_File_List ]

#------------------------------------------------------------------------------
# Define Date and Time Collector
#------------------------------------------------------------------------------

def Date_Time_Collector():
    
    Date = datetime.now().strftime("%Y_%m_%d")

    Time = datetime.now().strftime("%H_%M_%S")

    return Date, Time

#%%===========================================================================2
# 2. Helper Functions
#==============================================================================

class NpEncoder(json.JSONEncoder):
    """https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable"""
    def default(self, obj):
        if isinstance(obj, integer):
            return int(obj)
        if isinstance(obj, floating):
            return float(obj)
        if isinstance(obj, ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
def Energy_Wavelength_Converter( Energy_or_Wavelength ):

    """Convert from photon energy (in eV) to wavelength (in nm) or vice versa."""

    return h * c / ( e * Energy_or_Wavelength * 1e-9 )

def Radian_to_Degree_Converter( Angle_in_Radians ):
    
    """Convert an angle, or an array of angles, from units of radians to 
    
    units of degrees. This approach is around 40x quicker than importing and 
    
    using math.degrees"""
    
    return 180 * Angle_in_Radians / pi

def Degree_to_Radian_Converter( Angle_in_Degrees ):
    
    """Convert an angle, or an array of angles, from units of degrees to 
    
    units of radians. This approach is around 40x quicker than importing and 
    
    using math.radians"""
    
    return pi * Angle_in_Degrees / 40

def Photon_Flux_Spectrum_Loader( Spectrum_Name,

                                Light_Source_Path = Light_Source_Path ):
    
    """Using a spectrum path, load the spectrum's wavelengths and fluxes."""

    # --------------------------------------------------------------------------
    # Write Path to File
    # --------------------------------------------------------------------------

    Spectrum_Name = Spectrum_Name.replace(' ', '_')

    Spectrum_Path = path.join( Light_Source_Path, Spectrum_Name + '.xlsx')

    Spectrum_DataFrame = read_excel( Spectrum_Path )

    # --------------------------------------------------------------------------
    # Seperate Wavelengths and Fluxes
    # --------------------------------------------------------------------------

    Wavelengths = array( Spectrum_DataFrame[ Spectrum_DataFrame.columns[ 0 ] ] ) # nm
    
    Fluxes = array( Spectrum_DataFrame[ Spectrum_DataFrame.columns[ 1 ] ] ) # Photons/m2/s/nm 

    return Wavelengths, Fluxes

def Light_Power_Calculator( Wavelengths,

                           Energies,

                           Fluxes):
    
    """Determine the integrated light power (irradiance) of a spectrum using 
    
    its fluxes (photons/m2/s/nm) and the corresponding energies (in eV) and 
    
    wavelengths (in nm)."""

    # Units of Energies should be in eV, then multiply by e

    return simpson( e * Energies * Fluxes, x = Wavelengths ) / 10 # Factor of 10 converts from W/m2 to mW/cm2

def Planck_Photon_Flux_Wavelength( Wavelength, kT ):
    
    """At a given photon wavelength (in nm) and temperature (in K), determine 

    the photon flux as defined by the Planck spectrum in the limit that the 

    Maxwell-Boltzmann approximation is valid. Give this flux in units of 

    photons/m2/s/nm."""

    Pre_Factor = 2 * pi * c / ( 1e-9 * Wavelength ) ** 4 * 1e-9 # units are 1/s/m2/nm

    return Pre_Factor / ( exp( Energy_Wavelength_Converter( Wavelength ) / kT ) - 1 ) # Photons/m2/s/nm

def Lux_Constant_Calculator( Light_Source, Wavelengths, Energies ):
    
    """Determine the lux constant for a given light source, loading in the
    
    necessary luminous efficacy data."""
    
    #--------------------------------------------------------------------------
    # Load and normalise the photon flux spectrum
    #--------------------------------------------------------------------------

    Photon_Fluxes = Photon_Flux_Interpolator_and_Scaler(Light_Source,
                                                                 Wavelengths,
                                                                 1 ) # photons/m2/s/nm
        
    Photon_Irradiances = e * Energies * Photon_Fluxes # W/m2/nm
    
    P_light = simpson( Photon_Irradiances, x = Wavelengths )
        
    Normalised_Irradiances = Photon_Irradiances / P_light
    
    #--------------------------------------------------------------------------
    # Load interpolate the 2-deg luminous efficacy data
    #--------------------------------------------------------------------------
    
    V_Path = path.join( Supporting_Files_Path, "Illuminance_Files" )
    
    # V_Path = path.join( V_Path, "Luminous_Efficiency_Data.xlsx" )
    
    # V_DF = read_excel( V_Path, sheet_name = "2_Deg_V" )
    
    # Data from CIE: https://cie.co.at/datatable/cie-spectral-luminous-efficiency-photopic-vision
    
    V_Path = path.join( V_Path, "CIE_sle_photopic.csv" )

    V_DF = read_csv( V_Path )

    Luminous_Efficiencies = interp( Wavelengths, V_DF[ V_DF.columns[ 0 ] ], V_DF[ V_DF.columns[ 1 ] ] )
    
    return 683 * simpson( Luminous_Efficiencies * Normalised_Irradiances, x = Wavelengths ) # Constant of prop between lux and W/m2 in units of lux * W/m2
    
def Photon_Flux_Interpolator_and_Scaler( Spectrum_Name,

                                        New_Wavelengths,

                                        New_Light_Power,

                                        Light_Source_Path = Light_Source_Path):
    
    """Using a set of x-data and y-data, interpolate to determine what the 

    values would be at the desired x data points. Do this using the name of 

    the light (Light_Spectrum_Type) and the desired wavelengths 

    (New_Wavelengths) and the desired integrated intensity (New_Light_Power)

    in units of mW/cm2."""

    # --------------------------------------------------------------------------
    # Load in the spectral data and total irradiance
    # --------------------------------------------------------------------------

    Wavelengths, Fluxes = Photon_Flux_Spectrum_Loader(Spectrum_Name,

                                                      Light_Source_Path = Light_Source_Path) # nm, photons/m2/s/nm

    Energies = Energy_Wavelength_Converter( Wavelengths )  # in eV

    P_light = Light_Power_Calculator( Wavelengths,

                                     Energies,

                                     Fluxes) # mW/cm2

    # --------------------------------------------------------------------------
    # Determine bounds on the data
    # --------------------------------------------------------------------------

    min_wavelength = min(Wavelengths)

    max_wavelength = max(Wavelengths)

    # --------------------------------------------------------------------------
    # Define new photon fluxes as the new wavelengths
    # --------------------------------------------------------------------------

    New_Energies = Energy_Wavelength_Converter(New_Wavelengths)

    if P_light != 0:

        New_Fluxes = New_Light_Power * array([

            interp( wavelength,

                   Wavelengths,

                   Fluxes / P_light ) # Units are Units of fluxes (new_light_power and light_power have same units)

            if min_wavelength < wavelength < max_wavelength else 0

            for wavelength in New_Wavelengths])

    else:

        New_Fluxes = array( [ 0 for wavelength in New_Wavelengths ] )

    return New_Fluxes # photons/m2/s/nm

def EQE_Loader_and_Interpolator(Filename, Energies, Lower_Limit,

                                EQE_Spectra_Folder_Path=EQE_Spectra_Folder_Path):
    """Load and interpolate an EQE spectrum before trimming off the lower energies
    according to the thermodynamic lower limit (see https://www.nature.com/articles/s41467-020-19434-0
    and https://onlinelibrary.wiley.com/doi/full/10.1002/solr.202300277  for details)"""

    # --------------------------------------------------------------------------
    # First, Load Data
    # --------------------------------------------------------------------------

    Data = read_excel(path.join(EQE_Spectra_Folder_Path, Filename))

    Columns = Data.columns

    Wavelength_Data, EQE_Data = array(
        list(Data[Columns[0]])), array(list(Data[Columns[1]]))

    Energy_Data = Energy_Wavelength_Converter(Wavelength_Data)

    if max(EQE_Data) > 1:  # Convert from % to unitless
        EQE_Data /= 100

    # --------------------------------------------------------------------------
    # Next, Interpolate Data (Must Be in Terms of Increasing x)
    # --------------------------------------------------------------------------

    Reverse_Output = False

    if Energies[1] < Energies[0]:
        Reverse_Output = True
        Energies = Energies[::-1]

    if Energy_Data[1] < Energy_Data[0]:
        Energy_Data = Energy_Data[::-1]
        EQE_Data = EQE_Data[::-1]

    # --------------------------------------------------------------------------
    # Next, Interpolate Data (Must Be in Terms of Increasing x)
    # --------------------------------------------------------------------------

    EQEs = interp(Energies, Energy_Data, EQE_Data)

    for j, E in enumerate(Energies):
        if E < Lower_Limit:
            EQEs[j] = 0

    if Reverse_Output:
        EQEs = EQEs[::-1]

    EQEs = array( [ max( 0, j ) for j in EQEs ] ) # can't have -ve EQEs 
    
    return EQEs

def Optical_Constant_Loader_and_Interpolator( Wavelengths, Material_Type,
                                             
                                             Supporting_Files_Path = Supporting_Files_Path ):
    
    """For a given material type, import and interpolate its optical constants 
    
    at a given set of wavelengths."""
    
    Material_Path = path.join( Supporting_Files_Path, "Optical_Constants" )
    
#    print( Material_Type )
    
    Material_Path = path.join( Material_Path, Material_Type )
        
    Raw_ns = read_excel( Material_Path + '_n.xlsx' )  # first column is wavelengths, second is ns
    
    Raw_ks = read_excel( Material_Path + '_k.xlsx' )  # first column is wavelengths, second is ks
    
    ns = interp( Wavelengths, list( Raw_ns[ Raw_ns.columns[ 0 ] ] ), list( Raw_ns[ Raw_ns.columns[ 1 ] ] ) )
    
    ks = interp( Wavelengths, list( Raw_ks[ Raw_ks.columns[ 0 ] ] ), list( Raw_ks[ Raw_ks.columns[ 1 ] ] ) )
    
    return ns, ks

def Delta_V_oc_NR_Empirical( E_g,
                  
                  Delta_V_0 = 0.57, # V/eV
                  
                  eps = 1.49, # eV
                  
                  r = 10, # 1/eV
                  
                  delta = 0.16, # eV 
                  
                  ):
                      
    """Empirical model for non-radiative Voc loss in organic PVs from

    https://onlinelibrary.wiley.com/doi/10.1002/solr.202400456"""

    
    if E_g == eps:
        
        # treat special case first
        
        return Delta_V_0 * ( 1 / r + delta )
    
    else:
            
        denominator = 1 - exp( - r * ( eps - E_g ) )
        
        return Delta_V_0 * ( ( eps - E_g ) / denominator + delta )
                  
def Board_Efficiency_Calculator( V_mpp, I_mpp ):
    
    """Using the maximum power point voltage and current (not current density)
    
    determine the board efficiency using data for three current densities,
    
    then interpolate to find value at given current density."""
    
    #--------------------------------------------------------------------------
    # Interpolate in 0.1mA case:
    #--------------------------------------------------------------------------
    
    V_100uAs = array( I_100uA_Data[ I_100uA_Data.columns[0] ] )
    
    eta_100uAs = array( I_100uA_Data[ I_100uA_Data.columns[1] ] ) / 100
        
    if V_mpp < 0 or V_mpp > 3.6294 / 0.6660: # Root in linear extrapolation model
        
        eta_100uA = 0
        
    elif V_mpp < min( V_100uAs ):
        
        eta_100uA = 2.6665 * V_mpp + 0.1449 # Linear Extrpolation 
    
    elif V_mpp > max( V_100uAs ):
        
        eta_100uA = -0.6660 * V_mpp + 3.6294 # Linear Extrpolation 
        
        
    else: # Data within experimental range; interpolate:
        
        eta_100uA = interp( V_mpp, V_100uAs, eta_100uAs )
    
    #--------------------------------------------------------------------------
    # Interpolate in 1mA case:
    #--------------------------------------------------------------------------
    
    V_1mAs = array( I_1mA_Data[ I_1mA_Data.columns[0] ] )
    
    eta_1mAs = array( I_1mA_Data[ I_1mA_Data.columns[1] ] ) / 100
        
    if V_mpp < 0 or V_mpp > 3.7302 / 0.6845: # Root in linear extrapolation model
        
        eta_1mA = 0
        
    elif V_mpp < min( V_1mAs ):
        
        eta_1mA = 0.7861 * V_mpp + 0.5933
            
    elif V_mpp > max( V_1mAs ):
        
        eta_1mA = -0.6845 * V_mpp + 3.7302        
        
    else: # Data within experimental range; interpolate:
        
        eta_1mA = interp( V_mpp, V_1mAs, eta_1mAs )    
    
    #--------------------------------------------------------------------------
    # Interpolate in 10mA case:
    #--------------------------------------------------------------------------
    
    V_10mAs = array( I_10mA_Data[ I_10mA_Data.columns[0] ] )
    
    eta_10mAs = array( I_10mA_Data[ I_10mA_Data.columns[1] ] ) / 100
        
    if V_mpp < 0 or V_mpp > 4.4868 / 0.8634: # Root in linear extrapolation model
        
        eta_10mA = 0
        
    elif V_mpp < min( V_10mAs ):
        
        eta_10mA = 0.7415 * V_mpp + 0.5882
            
    elif V_mpp > max( V_10mAs ):
        
        eta_10mA = -0.8634 * V_mpp + 4.4868 
        
    else: # Data within experimental range; interpolate:
        
        eta_10mA = interp( V_mpp, V_10mAs, eta_10mAs )    
        
    #--------------------------------------------------------------------------
    # Interpolate BETWEEN Currents to Estimate eta at I_mpp
    #--------------------------------------------------------------------------
            
    ### Data is spread logarithmically 
    
    log_Is = [ log10( 0.1 ), log10( 1 ), log10( 10 ) ]
    
    etas = [ eta_100uA, eta_1mA, eta_10mA ]
    
    eta_interp = interp( log10( I_mpp ), log_Is, etas )
        
    return eta_100uA, eta_1mA, eta_10mA, eta_interp


#%%===========================================================================3
# 3. Functions for Small-Area Thermodynamic Limit Calculations
#------------------------------------------------------------------------------

def w(E, V, kT):
    """Determine the factor to correct for band-filling effects in a given 

    layer using the photon energy in eV (E), the voltage across that layer (V),

    and the thermal energy (kT) in eV."""
        
    return 1 / (1 + exp( (V - E) / 2 / kT)) ** 2

def Band_Filling_Factor(E, V, kT):
    """Determine the factor to correct EQE for band-filling effects in a given 

    layer using the photon energy in eV (E), the voltage across that layer (V),

    and the thermal energy (kT) in eV."""

    return tanh( (E - V) / 4 / kT )

def V_eff_Calculator( V, J, A, R_s ):
    
    """Determine effective voltage drop using actual voltage drop, V,
    
    current density (J in mA/cm2), Area (A in cm2), and series resistance
    
    (R_s, in Ohms)."""
    
    return V - J / 1000 * A * R_s # factor of 1000 to convert mA to A

def J_Calculator( V,
                   
                   EQEs, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL,
                   
                   kT,
                   
                   R_sh,
                   
                   R_s,
                   
                   A ):

    """For a given voltage V, evaluate the current density by calculating the 
    
    photocurrent density, the dark saturation current density, and the 
    
    contribution of the shunt term is optional (neglected if R_sh = inf).
    
    
    A non-zero series resistance is accounted for using V_eff and a recursive
    
    function."""
    
    #--------------------------------------------------------------------------
    # First, Treat Case when there is no series resistance:
    #--------------------------------------------------------------------------    
    
    if R_s == 0:
        
        #--------------------------------------------------------------------------
        # Calculate Photocurrent
        #--------------------------------------------------------------------------    
    
        Corrections_tanh = Band_Filling_Factor( Energies, V, kT )
        
        J_ph = e * simpson( EQEs * Corrections_tanh * Incident_Photon_Fluxes, x = Wavelengths ) # Units are C/m2/s = A/m2
    
        J_ph = J_ph / 10 # convert to mA/cm2
        
#        print( J_ph, e * simpson( EQEs[ ::-1] * Corrections_tanh[::-1] * Incident_Photon_Fluxes[::-1], x = Wavelengths ) )
        
        #--------------------------------------------------------------------------
        # Calculate Dark Saturation Current
        #--------------------------------------------------------------------------    
    
        Corrections_w = w( Energies, V, kT )
        
        J_0 = e * simpson( EQEs * Corrections_w * Black_Body_Fluxes, x = Wavelengths ) / EQE_EL # Units are A/m2
    
        J_0 = J_0 / 10 # convert to mA/cm2
        
        #--------------------------------------------------------------------------
        # Calculate Total Current Density
        #--------------------------------------------------------------------------    
    
        J = J_0 * ( exp( V / kT ) - 1 ) - J_ph
        
        #--------------------------------------------------------------------------
        # Add shunt Term if not Infinite
        #--------------------------------------------------------------------------    
        
        if R_sh != inf:
            
            J += V / R_sh * 1000 # Factor of 1000 to convert to mA/cm2
         
    #--------------------------------------------------------------------------
    # Otherwise, find J as root of expression.
    #--------------------------------------------------------------------------    
    
    else:    
        
        def Function_to_Find_Root_of( J ):
            
            # USE RECURSIVE CODING
            
            V_eff = V_eff_Calculator( V, J, A, R_s )
            
            J_eff = J_Calculator( V_eff,
                               
                               EQEs, 
                               
                               Incident_Photon_Fluxes,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL,
                               
                               kT,
                               
                               R_sh,
                               
                               0, # Let Rs = 0 as V_eff includes it
                               
                               A )  # returns J_eff in mA/cm2
            
            return J_eff - J # Goes to nought when solution found
        
        J_bound = 0.01
        
        Bound_Found = False
        
        while not Bound_Found:
            
            try:
                
                brentq( Function_to_Find_Root_of, -J_bound, J_bound )
                
            except:
                
                J_bound *= 10
                
            else:
                
                Bound_Found = True
            
        J = brentq( Function_to_Find_Root_of, -J_bound, J_bound )
        
    return J # units are mA/cm2

def V_oc_Finder( EQEs, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL,
                   
                   kT,
                   
                   R_sh ):
    
    """Find the V_oc across a cell using the incident photon fluxes,
    
    EQEs, R_sh, etc. The series resistance becomes negligible."""
    
    def Function_to_Find_Root_of( V ):
        
        return J_Calculator( V,
                           
                           EQEs, 
                           
                           Incident_Photon_Fluxes,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL,
                           
                           kT,
                           
                           R_sh,
                           
                           0,
                           
                           1 ) # units are mA/cm2
    
    Max_V = 1
    Bound_Found = False
    
    while not Bound_Found:
        
        try:
            
            brentq( Function_to_Find_Root_of, 0, Max_V )
            
        except:
            
            Max_V += 0.2
            
        else:
            
            Bound_Found = True
                
    return brentq( Function_to_Find_Root_of, 0, Max_V )

def V_mpp_Calculator(
                   
                   EQEs, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL,
                   
                   kT,
                   
                   R_sh,
                   
                   R_s,
                   
                   A,
                   
                   Min_V = 0,
                   
                   Max_V = 1 ):
    
    def Function_to_Minimize( V ):
        
        return V * J_Calculator( V,
                           
                           EQEs, 
                           
                           Incident_Photon_Fluxes,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL,
                           
                           kT,
                           
                           R_sh,
                           
                           R_s,
                           
                           A ) # units are mW/cm2

    x0 = 0.8 * Max_V
                
    return minimize( Function_to_Minimize, x0, bounds = ( ( Min_V, Max_V ), ) )[ 'x' ][ 0 ] # units are mW/cm2

def Small_Area_Figure_of_Merit_Calculator( Light_Power,
        
                   EQEs, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL,
                   
                   kT,
                   
                   R_sh,
                   R_s,
                   
                   Area_cell ):
    
    """Evaluate the photovoltaic figures-of-merit for a small-area device,
    
    using the EQEof each cell, the incident photon fluxes (in photons/m2/s/nm), 
    
    the wavelengths (in nm), the light power/irradiance (in mW/cm2),
    
    and the thermal energy (kT, in eV)."""
            
    #--------------------------------------------------------------------------
    # Find J_sc
    #--------------------------------------------------------------------------    
           
    J_sc = abs( J_Calculator( 0,
                       
                       EQEs, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL,
                       
                       kT,
                       
                       R_sh,
                       
                       R_s,
                       
                       Area_cell ) ) # units are mA/cm2
    
    #--------------------------------------------------------------------------
    # Find V_oc
    #--------------------------------------------------------------------------    
            
    V_oc = V_oc_Finder( EQEs, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL,
                       
                       kT,
                       
                       R_sh ) # units are V
 
   # print( R_sh, V_oc )
    
    #--------------------------------------------------------------------------
    # Find V_mpp
    #--------------------------------------------------------------------------    

    V_mpp = V_mpp_Calculator(
                       
                       EQEs, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL,
                       
                       kT,
                       
                       R_sh,
                       
                       R_s,
                       
                       Area_cell,
                       
                       Min_V = 0,
                       
                       Max_V = V_oc ) # units are V
    
    #--------------------------------------------------------------------------
    # Find J_mpp
    #--------------------------------------------------------------------------    
    
    J_mpp = J_Calculator( V_mpp,
                       
                       EQEs, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL,
                       
                       kT,
                       
                       R_sh,
                       
                       R_s,
                       
                       Area_cell ) # units are mA/cm2
    
    #--------------------------------------------------------------------------
    # Calculate and Return Remaining Figures-of-Merit
    #--------------------------------------------------------------------------    
        
    P_mpp = abs( V_mpp * J_mpp ) # mW/cm2
    
    PCE = P_mpp / Light_Power
        
    return { "J_sc" : J_sc,
            
            "V_oc" : V_oc,
            
            "V_mpp" : V_mpp,
            
             "J_mpp" : J_mpp,
             
             "P_mpp" : P_mpp,
             
             "PCE" : PCE }

#%%===========================================================================4
# 4. Functions for Simulating EQE Spectra
#------------------------------------------------------------------------------

def Step_Function_EQE_Simulator( Energies, Energetic_Gap, B, A ):
    
    """Simulate a step function EQE spectrum using a list of energies, a value 
    
    for the above-gap EQE (A), a value for the below-gap EQE (B). In the SQ model,
    
    B = 0 and A = 1"""
                
    return array( [ A if E > Energetic_Gap else B for E in Energies ] )

def Urbach_Tail_EQE_Simulator( Energies, Energetic_Gap, Urbach_Energy, A ):
    
    """Simulate an EQE spectrum for a list of energies (in eV) using values for the 
    
    above-gap EQE (A) and the Urbach energy (in eV)."""
    
    return A * array( [ 1 if E > Energetic_Gap else 
                                   
                                   exp( ( E - Energetic_Gap ) / Urbach_Energy ) for E in Energies ] )


def Exction_Absorption_EQE_Simulator( Energies, Energetic_Gap, sigma_s, A, kT ):
    
    """Simulate an EQE spectrum for a list of energies (in eV) using values for the 
    
    above-gap EQE (A), the energetic gap (in eV), the energetic disorder (sigma_s, 
                                                                          
    in eV), and the thermal energy, kT (k * temperature / e, in units of eV)."""
        
    if sigma_s != 0:
        
        Exponential_Term = exp( ( Energies - Energetic_Gap + sigma_s ** 2 / 2 / kT  ) / kT )

        Erfc_Term = erfc( ( Energies - Energetic_Gap + sigma_s ** 2 / kT  ) / sigma_s / 2 ** 0.5 )
        
        Erf_Term_1 = erf( ( Energetic_Gap ) / sigma_s / 2 ** 0.5 )
        
        Erf_Term_2 = erf( ( Energies - Energetic_Gap ) / sigma_s / 2 ** 0.5 ) 
        
        EQEs = A / 2 * ( Exponential_Term * Erfc_Term + Erf_Term_1 + Erf_Term_2 )
        
    else:
        
        EQEs = Urbach_Tail_EQE_Simulator( Energies, Energetic_Gap, kT, A )
        
    return EQEs


#%%===========================================================================5
# 5. Optical Modelling
#==============================================================================

"""Summary:

In this section functions are defined for determining the Fresnel reflection/
transmission coefficients for all layers in the device architecture (using
the respective optical constants for each layer of the device and the 
surrounding medium). 

Following this, the optical electric field distribution is determined through
optical transfer matrix modelling.

In addition, the reflectivity of any encapsulating later and/or the albedo of 
any surface on the posterior side of the photovoltaic is accounted for using
the approach outlined in the supporting information of:
    
https://onlinelibrary.wiley.com/doi/10.1002/solr.202400456

This section of the script expands on the work of Burkhard et al.:

G. F. Burkhard, E. T. Hoke, M. D. McGehee, Adv. Mater. 2010, 22, 3293.
"""

def phi_j( Refractive_Index_0 , Refractive_Index_j , Incident_Angle ):
    
    """Determine the angle of propagation in layer j in radians using the 
    
    incident angle in layer 0 (also in radians), the refractive index in the 
    
    zeorth layer, and the refractive index in the j-th layer."""
    
    return arcsin( Refractive_Index_0 / Refractive_Index_j * sin( Incident_Angle ) )

def q_j_calculator( Complex_Refractive_Index, Angle_of_Propagation ):
    
    """Determine the quantity q_j for the j-th layer, using its complex 
    
    refractive index and the angle of propagation (in radians)."""
    
    return Complex_Refractive_Index * cos( Angle_of_Propagation )

def r_jk_s( q_j , q_k ):
    
    """Determine the Fresnel reflection coeffient for s-polarised light."""
    
    return ( q_j - q_k ) / ( q_j + q_k )

def t_jk_s( q_j , q_k ):
    
    """Determine the Fresnel transmission coeffient for s-polarised light."""
    
    return 2 * q_j / ( q_j + q_k )

def r_jk_p( q_j , q_k , Complex_Refractive_Index_j , Complex_Refractive_Index_k ):
    
    """Determine the Fresnel reflection coeffient for p-polarised light."""
    
    Numerator = Complex_Refractive_Index_j ** 2 * q_k - Complex_Refractive_Index_k ** 2 * q_j 
    
    Denominator =  Complex_Refractive_Index_k ** 2 * q_j + Complex_Refractive_Index_j ** 2 * q_k
    
    return Numerator / Denominator

def t_jk_p( q_j , q_k , Complex_Refractive_Index_j , Complex_Refractive_Index_k ):
    
    """Determine the Fresnel reflection coeffient for p-polarised light."""
    
    Numerator =  2 * Complex_Refractive_Index_j * Complex_Refractive_Index_k * q_j
    
    Denominator =  Complex_Refractive_Index_k ** 2 * q_j + Complex_Refractive_Index_j ** 2 * q_k
    
    return Numerator / Denominator

def Attenuation_Delta( alpha, d, phi ):
    
    """Calculate the intensity attenutation during a pass through the substrate 
    
    with absorption coefficient alpha and thickness d (the units of these 
                                                       
    parameters must be inverse). The angle of propagation in the substrate is 
    
    expected to be previously calculated, and it should be entered in radians."""
    
    return exp( - alpha * d * cos( phi ) )

def R_Intensity( r_Fresnel ):
    
    """Determine the (real) reflection coefficient for the incident light 
    
    intensity, using the complex Fresnel reflection coefficient."""
    
    return abs( r_Fresnel ) ** 2

def T_Intensity( t_Fresnel, eta_j, eta_k ):

    """Determine the (real) transmission coefficient for the incident light
    
    intensity, using the complex Fresnel transmission coefficient and the 
    
    (real) refractive indices of the two media (with eta_k being the refractive
                                                
    index of the final medium)."""
    
    return eta_k / eta_j * abs( t_Fresnel ) ** 2

def T_trans( R_01, T_01, R_stack, Delta_1 ):
    
    """Determine the intensity transmission of light up to the interface 
    
    between the incoherent substrate and the thin-film stack, using the 
    
    reflection and transmission coefficients of the first interface (between 
                                                                     
    the surrounding medium and the substrate), the reflectance of the thin-film 
    
    stack (determined through optical modelling techniques), and the 
    
    attenuation in the substrate (Delta_1)."""
    
    return T_01 * Delta_1 / ( 1 - Delta_1 ** 2 * R_01 * R_stack )    

def R_Tot( R_01, R_stack, Delta_1 ):
    
    """Determine the fraction of light reflected from the device using the 
    
    reflection coefficients of the first interface (between the surrounding 
                                                    
    medium and the substrate), the reflectance of the thin-film stack 
    
    (determined through optical modelling techniques), and the attenuation in 
    
    the substrate (Delta_1)."""
    
    return ( R_01 + Delta_1 ** 2 * ( 1 - 2 * R_01 ) * R_stack ) / ( 1 - Delta_1 ** 2 * R_01 * R_stack )    
            
def T_Tot( R_01, R_stack, Delta_1 ):
    
    """Determine the fraction of light transmitted through the device, using 
    
    the reflection coefficients of the first interface (between the surrounding 
                                                        
    medium and the substrate), the reflectance of the thin-film stack (determined 
    
    through optical modelling techniques), and the attenuation in the substrate 
                                                                       
    (Delta_1)."""
    
    return ( ( 1 - R_01 ) * Delta_1 * ( 1 - R_stack ) ) / ( 1 - Delta_1 ** 2 * R_01 * R_stack )    

def r_Fresnel( Transfer_Matrix ):
    
    """Using the 2x2 transfer matrix, determine the complex Fresnel reflection
    
    coefficient as the ratio of the 21 element to the 11 element. This can be 
    
    done for the total transfer matrix and the partial transfer matrices."""
    
    return Transfer_Matrix[ 1, 0 ] / Transfer_Matrix[ 0, 0 ]

def t_Fresnel( Transfer_Matrix ):
    
    """Using the 2x2 transfer matrix, determine the complex Fresnel 
    
    transmission coefficient as the inverse of 11 element. This can be done for 
    
    the total transfer matrix and the partial transfer matrices."""
    
    return 1 / Transfer_Matrix[ 0, 0 ]

def I_mat( n_tilde1, n_tilde2 ):
    
    """Compute the interface matrix using the complex refractive indices of the 
    
    different layers (labelled as n1 and n2)."""
    
    r = ( n_tilde1 - n_tilde2 ) / ( n_tilde1 + n_tilde2 )
    
    t = ( 2 * n_tilde1 ) / ( n_tilde1 + n_tilde2 )
    
    return array( [ [ 1, r ],
                    [ r, 1 ] ] ) / t

def L_mat( n_tilde, d, Wavelength ):
    
    """Determine the propgation matrix for a given layer, using its complex
    
    refractive index (n_tilde), thickness (d), and the wavelength of the light."""
    
    xi = 2 * pi * d * n_tilde / Wavelength
    
    return array( [ [ exp( - 1j * xi ), 0 ], 
                   
                    [ 0, exp( 1j * xi ) ] ] )

def Single_Pass_Optical_Modeller( Energies, Wavelengths, IQEs, Photon_Fluxes, 
                                 
    Surrounding_Medium, Layer_Types, Layer_Thicknesses, Position_Spacing,
                     
    Active_Layer_Index,
                     
    Incoherent_Substrate = True,
                    
    Infinite_Final_Layer = False,
    
    Substrate_Index = 1 ):
    

    """Conduct the optical modelling simulations using the wavelengths, the 
    
    internal quantum efficiencies (IQEs), the surrounding medium type, the 
    
    layer types and thicknesses, the separation between points, and the index
    
    of the active layer. The wavelengths are assumed to be in an list/array in 
    
    terms of increasing wavelength."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # If Final Layer is Finite, Add Surrounding Medium to Device
    #----------------------------------------------------------------------------------------------------------------------- 
    
    if not Infinite_Final_Layer:
    
        Layer_Types.append( Surrounding_Medium )
    
        Layer_Thicknesses.append( 0 )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Optical Constants
    #----------------------------------------------------------------------------------------------------------------------- 
    
    n_tildes = zeros( ( len( Layer_Types ), len( Wavelengths ) ), dtype = complex )

    # load index of refraction for each material in the stack

    for i, Layer_Type in enumerate( Layer_Types ):
        
#        print( i, Layer_Type, Layer_Types )
        
        n_is, k_is = Optical_Constant_Loader_and_Interpolator( Wavelengths, Layer_Type )

        n_tildes[ i, : ] = n_is + 1j * k_is

    #----------------------------------------------------------------------------------------------------------------------
    # Determine Incoherent Intensity Transmission Through Substrate
    #----------------------------------------------------------------------------------------------------------------------

    n_0s, k_0s = Optical_Constant_Loader_and_Interpolator( Wavelengths, Surrounding_Medium ) 

    n_0_tildes = n_0s + 1j * k_0s  

    n_1_tildes = n_tildes[ 0, : ]

    r_01s = ( n_0_tildes - n_1_tildes ) / ( n_0_tildes + n_1_tildes )
    
    t_01s = 2 * n_0_tildes / ( n_0_tildes + n_1_tildes )

    R_01s = abs( r_01s ) ** 2

    T_01s = abs( n_1_tildes / n_0_tildes ) * abs( t_01s ) ** 2
    
    #----------------------------------------------------------------------------------------------------------------------
    # Simulate Positions
    #----------------------------------------------------------------------------------------------------------------------

    if Incoherent_Substrate:
        
        Substrate_Thickness = Layer_Thicknesses[ 0 ]
    
        Layer_Thicknesses[ 0 ] = 0
    
    Positions = []

    for j, Type in enumerate( Layer_Types ):
    
        Positions_j = arange( 0, Layer_Thicknesses[ j ], Position_Spacing )[ 1: ]
    
        if len( Positions_j ) > 0:
        
            Positions.append( Positions_j )

        else:
        
            Positions.append( array( [ 0 ] ) )
                
    #----------------------------------------------------------------------------------------------------------------------
    # Calculate Transfer Matrices, and Normalised Optical Electric Field at Each Wavelength and Position
    #----------------------------------------------------------------------------------------------------------------------

    # Intialise Rs and Ts, and Absorptances
    
    global R_tots, T_tots

    R_stacks = zeros( len( Wavelengths ) )  # Reflectance of Thin-Film Stack

    T_stacks = zeros( len( Wavelengths ) )  # Transmittance of Thin-Film Stack

    R_tots = zeros( len( Wavelengths ) )  # Total Reflectance of Device
    
    T_tots = zeros( len( Wavelengths ) )  # Total Transmittance of Device
    
    T_subs = zeros( len( Wavelengths ) )
    
    if Incoherent_Substrate:
                
        Deltas = zeros( len( Wavelengths ) )
        
    E2s = [ zeros( ( len( Positions_j ), len( Wavelengths ) ) ) for Positions_j in Positions ]
                    
    # Iterate Over All Wavelengths
    
    Lower_Index = 1
    
    if not Incoherent_Substrate:
        
        Lower_Index = 0        

    for j, Wavelength in enumerate( Wavelengths ):
    
        #------------------------------------------------------------------------------------------------------------------
        # Calculate the Transfer Matrices for Reflection/Transmission of Stack
        #------------------------------------------------------------------------------------------------------------------

        Ordered_Matrices = [ I_mat( n_tildes[ Layer_Index, j ], n_tildes[ Layer_Index + 1, j ] ) 
                          
                          for Layer_Index in range( len( Layer_Thicknesses )  - 1 ) ]

        Shift = 0
    
        for Layer_Index in range( 1, len( Layer_Thicknesses ) ):
        
            Layer_Thickness = Layer_Thicknesses[ Layer_Index ]
        
            n_tilde_j = n_tildes[ Layer_Index, j ]
                
            Ordered_Matrices.insert( Layer_Index + Shift, L_mat( n_tilde_j, Layer_Thickness, Wavelength ) )
    
            Shift += 1
    
        # If the first layer is NOT incoherent, add the Interface matrix and layer matrix for the surrounding_medium
    
        Layer_1_Thickness = Substrate_Thickness
    
        if not Incoherent_Substrate:
                
            n_0_tilde = n_0_tildes[ j ]  

            n_1_tilde = n_tildes[ 0, j ]
        
            I = I_mat( n_0_tilde, n_1_tilde )
                
            L = L_mat( n_1_tilde, Layer_1_Thickness, Wavelength )
    
            Ordered_Matrices = [ I, L ] + Ordered_Matrices
    
        # Determine TOTAL System Transfer Matrix
        
        S = linalg.multi_dot( Ordered_Matrices )
    
        # Determine Reflection and Transmission of Thin-Film Stack

        R_stack = abs( S[ 1, 0 ] / S[ 0, 0 ] ) ** 2
        
        R_stacks[ j ] = R_stack
        
        if Incoherent_Substrate:
            
            T_stack = abs( n_0_tildes[ j ] / n_tildes[ 0, j ] ) * abs( 1 / S[ 0, 0 ] ) ** 2
            
            T_stacks[ j ] = T_stack 
            
            Delta = exp( - 4 * pi * n_tildes[ 0, j ].imag / Wavelength * Layer_1_Thickness )
            
            Deltas[ j ] = Delta
                        
            T_subs[ j ] = T_01s[ j ] * Delta / ( 1 - Delta ** 2 * R_01s[ j ] * R_stacks[ j ] )      
            
            R_tots[ j ] = R_01s[ j ] + T_01s[ j ] ** 2 * Delta ** 2 * R_stacks[ j ] / ( 1 - Delta ** 2 * R_01s[ j ] * R_stacks[ j ] )
    
            T_tots[ j ] = T_subs[ j ] * T_stack
        
        else:
            
            T_stack = abs( 1 / S[ 0, 0 ] ) ** 2
            
            T_stacks[ j ] = T_stack
            
            R_tots[ j ] = R_stack
            
            T_tots[ j ] = T_stack
        
            T_subs[ j ] = 1
        
        #------------------------------------------------------------------------------------------------------------------
        # Calculate Anterior and Posterior Partial System Transfer Matrices 
        #------------------------------------------------------------------------------------------------------------------

        for Layer_Index in range( Lower_Index, len( Layer_Thicknesses ) - 1 ):

            xi = 2 * pi * n_tildes[ Layer_Index, j ] / Wavelength

            Layer_Thickness = Layer_Thicknesses[ Layer_Index ]

            x = Positions[ Layer_Index ]      
        
            #--------------------------------------------------------------------------------------------------------------
            # Calculate S_Prime
            #--------------------------------------------------------------------------------------------------------------        
        
            if Incoherent_Substrate:
        
                Anterior_Matrices = Ordered_Matrices[ : 2 * Layer_Index - 1 ]
            
            else:
            
                Anterior_Matrices = Ordered_Matrices[ : 2 * Layer_Index + 1 ]
            
            if len( Anterior_Matrices ) > 1:

                S_prime = linalg.multi_dot( Anterior_Matrices )
                        
            else:
            
                S_prime = Anterior_Matrices[ 0 ]

            #--------------------------------------------------------------------------------------------------------------
            # Calculate S_pprime (double prime)
            #--------------------------------------------------------------------------------------------------------------                            
        
            if Incoherent_Substrate:
            
                Posterior_Matrices = Ordered_Matrices[ 2 * Layer_Index : ]
            
            else:
            
                Posterior_Matrices = Ordered_Matrices[ 2 * Layer_Index + 2 : ]
            
            if len( Posterior_Matrices ) > 1:
            
                S_pprime = linalg.multi_dot( Posterior_Matrices )
            
            else:
            
                S_pprime = Posterior_Matrices[ 0 ]
                
            num = sqrt( T_subs[ j ] ) * ( S_pprime[ 0, 0 ] * exp( -1j * xi *( Layer_Thickness - x ) ) + 
                          
                          S_pprime[ 1, 0 ] * exp( 1j * xi * ( Layer_Thickness - x ) ) )

            den = S_prime[ 0, 0 ] * S_pprime[ 0, 0 ] * exp( -1j * xi * Layer_Thickness ) + S_prime[ 0, 1 ] * S_pprime[ 1, 0 ] * exp( 1j * xi * Layer_Thickness )
            
            E2s[ Layer_Index ][ :, j ] = abs( n_0_tildes[ j ] / n_1_tildes[ j ] ) * abs( num / den ) ** 2
            
    #----------------------------------------------------------------------------------------------------------------------
    # Determine Absorption
    #----------------------------------------------------------------------------------------------------------------------
            
    Abs = {}
    
    for Layer_Index in arange( Lower_Index, len( Layer_Thicknesses ) ):
            
        alphas = 4 * pi * n_tildes[ Layer_Index, : ].imag / Wavelengths
        
        xs = Positions[ Layer_Index ]
        
        Norm_E_sqrs = E2s[ Layer_Index ]
        
        Abs_j = []
                            
        for j, Wavelength in enumerate( Wavelengths ):
            
            Es = Norm_E_sqrs[ :, j ]
            
            AbsRate =  alphas[ j ] * n_tildes[ Layer_Index, j ].real * Es
            
            Abs_j.append( sum( AbsRate ) * Position_Spacing )

        Abs[ 'Abs_' + str( Layer_Index + 1 ) + '_' + Layer_Types[ Layer_Index ] ] = Abs_j

    # Determine total absorption
    
    Abs_tot = 1 - R_tots - T_tots
    
    #----------------------------------------------------------------------------------------------------------------------
    # Determine EQE
    #----------------------------------------------------------------------------------------------------------------------
    
    EQEs = array( Abs[  'Abs_' + str( Active_Layer_Index ) + '_' + Layer_Types[ Active_Layer_Index - 1 ] ] ) * array( IQEs )
    
    #----------------------------------------------------------------------------------------------------------------------
    # Determine Generation Rate
    #----------------------------------------------------------------------------------------------------------------------
        
    # First, determine irradiances:
    
    Irradiances = Photon_Fluxes * Energies * e    
        
    # Next, determine energy dissipation
    
    alpha_ALs = 4 * pi * n_tildes[ Active_Layer_Index, : ].imag / Wavelengths
    
    Q_js_init = alpha_ALs * n_tildes[ Active_Layer_Index, : ].real / n_0_tildes.real * Irradiances   # Need to multiply by |E|^2

    # Calculate G_light
        
    G_lights = []
    
    E2s_Al = E2s[ Active_Layer_Index ]
    
    # Determine Wavelength Separation and Positions
    
    Delta_Wavelength = Wavelengths[ 1 ] - Wavelengths[ 0 ] 
        
    Active_Layer_Positions = Positions[ Active_Layer_Index ]

    for k, Position in enumerate( Active_Layer_Positions ):
                           
        Q_js = E2s_Al[ k, : ] * Q_js_init
                                
        G_lights_j = IQEs * Q_js / h / c * ( Wavelengths * 1e-9 ) 
            
        G_light = sum( G_lights_j ) * Delta_Wavelength # Units are excitons/s/cm2/nm
        
        G_light *= 1e7 # Units are excitons/s/cm3
                
        G_lights.append( G_light )
        
    #----------------------------------------------------------------------------------------------------------------------
    # Output Dictionaries
    #---------------------------------------------------------------------------------------------------------------------
    
    # Optical Electric Field:
    
    Electric_Field_Output = { 'Positions' : Positions, 'E2s' : E2s }
    
    # Spectral Parameters:
    
    Optical_Modelling_Output = { 'Wavelength' : Wavelengths,
                                
                                 'Energy' : Energies,
                                
                                 'EQE' : EQEs,
                                
                                 'R' : R_tots,
                                
                                 'T' : T_tots,
                               
                                 'Total Abs' : Abs_tot }
    
    Optical_Modelling_Output = Optical_Modelling_Output | Abs
    
    # Generation Rate:
        
    Generation_Rate_Output = { 'x' : Active_Layer_Positions,
                             
                               'G_light' : G_lights }        

    return Electric_Field_Output, Optical_Modelling_Output, Generation_Rate_Output

#%%===========================================================================6
# 6. Multi-Junction Simulation
#==============================================================================

def MJ_J_i_Calculator( V_i_prime,
                   
                   EQEs_i, 
                   
                   Incident_Photon_Fluxes_i,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_i,
                   
                   kT,
                   
                   R_sh_i,
                                      
                   Coupling = False,
                   
                   EQEs_em = 1,
                   
                   V_em = 0 ):

    """For the i-th cell in the mult-junction stack, evaluate the photocurrent
    
    density, the dark saturation current density, any photocurrent attributed
    
    due to radiative coupling, and ultimately the total photocurrent.
    
    
    The contribution of the shunt term is optional (neglected if R_sh_i = inf).
    
    In addition, the effect of radiative coupling can be switched off by setting
    
    'Coupling' = False. To calculate the coupling terms, the EQEs, emission 
    
    voltage (V_em) and electroluminescent external quantum efficiency of the
    
    emitting cell is also required."""
    
    #--------------------------------------------------------------------------
    # Calculate Photocurrent
    #--------------------------------------------------------------------------    

    Corrections_tanh = Band_Filling_Factor( Energies, V_i_prime, kT )
    
    J_ph_i = e * simpson( EQEs_i * Corrections_tanh * Incident_Photon_Fluxes_i, x = Wavelengths )

    J_ph_i = J_ph_i / 10 # convert to mA/cm2
    
    #--------------------------------------------------------------------------
    # Calculate Dark Saturation Current
    #--------------------------------------------------------------------------    

    Corrections_w = w( Energies, V_i_prime, kT )
    
    J_0_i = e * simpson( EQEs_i * Corrections_w * Black_Body_Fluxes, x = Wavelengths ) / EQE_EL_i

    J_0_i = J_0_i / 10 # convert to mA/cm2
    
    #--------------------------------------------------------------------------
    # Calculate Total Current Density
    #--------------------------------------------------------------------------    

    J = J_0_i * ( exp( V_i_prime / kT ) - 1 ) - J_ph_i
    
    #--------------------------------------------------------------------------
    # Add shunt Term if not Infinite
    #--------------------------------------------------------------------------    
    
    if R_sh_i != inf:
        
        J += V_i_prime / R_sh_i * 1000 # Factor of 1000 to convert to mA/cm2
        
    #--------------------------------------------------------------------------
    # Add Radiative Coupling Term if True
    #--------------------------------------------------------------------------    
            
    if Coupling:
        
        # Corrections_tanh already evaluated above for absorbing cell
                
        Corrections_w_emitting_cell = w( Energies, V_em, kT )
        
        J_0_12 = e * simpson( EQEs_i * Corrections_tanh * EQEs_em * Corrections_w_emitting_cell * Black_Body_Fluxes, x = Wavelengths )
        
        J_0_12 = J_0_12 / 10 # convert to mA/cm2 
        
        J -= J_0_12 / 2 * ( exp( V_em / kT ) - 1 )
            
    return J

def MJ_V_1_Calculator( V,
                   
                   EQEs_1, 
                   EQEs_2, 

                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                   
                   R_s_1,
                   R_s_2,
                   Area,
                   
                   Iteration = 0,
                   
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                                    
                   V_em_1 = 0.5,
                   V_em_2 = 0.5,
                   
                   ):
 
    AR_s_tot = Area * ( R_s_1 + R_s_2 )    
    
    def Function_to_Find_Root_of( V_1_prime ):
        
        #======================================================================
        # If First Iteration        
        #======================================================================        
        
        V_1_em_guess = 0
        V_2_em_guess = 0
        
        if Iteration == 0:
            
            #------------------------------------------------------------------
            # Evaluate Current Generated by the Top Cell
            #------------------------------------------------------------------
            
            J_1 = MJ_J_i_Calculator( V_1_prime,
                               
                               EQEs_1, 
                               
                               Incident_Photon_Fluxes.copy(),
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_1,
                               
                               kT,
                               
                               R_sh_1,
                                                  
                               Coupling = Coupling_2_to_1,
                               
                               EQEs_em = EQEs_2,
                               
                               V_em = V_2_em_guess  
                               
                               )
            
            #------------------------------------------------------------------
            # Evaluate Voltage on the Top Cell
            #------------------------------------------------------------------
            
            V_2_prime = V - V_1_prime - ( J_1 * AR_s_tot / 1000 ) # Factor of 1000 to conver from mA Ohm to V
                    
            #------------------------------------------------------------------
            # Evaluate Fluxes Indicent on the Bottom Cell
            #------------------------------------------------------------------
            
            Incident_Photon_Fluxes_2 = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_1_prime, kT ) ) * Incident_Photon_Fluxes.copy()
            
            #------------------------------------------------------------------
            # Evaluate Current Generated by the Bottom Cell
            #------------------------------------------------------------------
            
            J_2 = MJ_J_i_Calculator( V_2_prime,
                               
                               EQEs_2, 
                               
                               Incident_Photon_Fluxes_2,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_2,
                               
                               kT,
                               
                               R_sh_2,
                                                  
                               Coupling = Coupling_1_to_2,
                               
                               EQEs_em = EQEs_1,
                               
                               V_em = V_1_em_guess 
                               
                               )
            
        #======================================================================
        # If Not First Iteration:
        #======================================================================        
                    
        else:
            
            #------------------------------------------------------------------
            # Evaluate Current Generated by the Top Cell
            #------------------------------------------------------------------

            J_1 = MJ_J_i_Calculator( V_1_prime,
                               
                               EQEs_1, 
                               
                               Incident_Photon_Fluxes.copy(),
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_1,
                               
                               kT,
                               
                               R_sh_1,
                                                  
                               Coupling = Coupling_2_to_1,
                               
                               EQEs_em = EQEs_2,
                               
                               V_em = V_em_2  
                               
                               )
            
            #------------------------------------------------------------------
            # Evaluate Voltage on the Top Cell
            #------------------------------------------------------------------
            
            V_2_prime = V - V_1_prime - ( J_1 * AR_s_tot / 1000 ) # Factor of 1000 to conver from mA Ohm to V
                    
            #------------------------------------------------------------------
            # Evaluate Fluxes Indicent on the Bottom Cell
            #------------------------------------------------------------------
            
            Incident_Photon_Fluxes_2 = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_1_prime, kT ) ) * Incident_Photon_Fluxes.copy()
            
            #------------------------------------------------------------------
            # Evaluate Current Generated by the Bottom Cell
            #------------------------------------------------------------------
            
            J_2 = MJ_J_i_Calculator( V_2_prime,
                               
                               EQEs_2, 
                               
                               Incident_Photon_Fluxes_2,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_2,
                               
                               kT,
                               
                               R_sh_2,
                                                  
                               Coupling = Coupling_1_to_2,
                               
                               EQEs_em = EQEs_1,
                               
                               V_em = V_em_1 
                               
                               )
                 
        #----------------------------------------------------------------------
        # Return Difference
        #----------------------------------------------------------------------
        
        return J_1 - J_2
    
    #--------------------------------------------------------------------------
    # Find root of above function to determine V_1_prime
    #--------------------------------------------------------------------------
            
    return brentq( Function_to_Find_Root_of, 0, V )
        
def MJ_V_em_i_Calculator( J,
                                       
                       EQEs_i, 
                                              
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_i,
                       
                       kT,
                       
                       R_sh_i ):
    
    """Find the effective emission voltage of a cell under illumination by 
    
    treating it as an LED in the dark. Radiative coupling is neglected in this
    
    case."""
    
    # Emission voltage is INDEPENDENT of series resistance. It only depends on
    
    # the current flowing through the cell (J), not the voltage. 
    
    # Below solves for V_em_prime, which is the effective emission voltage of 
    
    # the cell (internally) regardless of any Ohmic loss.

    def Function_to_Find_Root_of( V_em_prime ):
        
        J_dark = MJ_J_i_Calculator( V_em_prime,
                           
                           EQEs_i, 
                           
                           0, # No incident fluxes
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL_i,
                           
                           kT,
                           
                           R_sh_i,
                           
                           Coupling = False,
                           
                           EQEs_em = 1,
                           
                           V_em = 0 )
        
        return J_dark - abs( J )
    
    Max_V = 1
    Bound_Found = False
    
    while not Bound_Found:
        
        try:
            
            brentq( Function_to_Find_Root_of, 0, Max_V )
            
        except:
            
            Max_V += 1
            
        else:
            
            Bound_Found = True
                        
    return brentq( Function_to_Find_Root_of, 0, Max_V )

def MJ_J_tot_Calculator( V,
                   
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                   
                   R_s_1,
                   R_s_2,
                   Area,
                  
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Convergence_Criterion = 1e-6,
                   Max_Iterations = 300 ): 

    """At a given voltage (V), evaluate the total current flowing through the 
    
    device using an iterative process in which the effective emission voltages
    
    of cell 1 and cell 2 are evaluated.
    
    
    The iterative process proceeds as follows:
        
    (1) a value for V_1 is identified by setting J_1(V_1) = J_2(V - V_1).
    
    (2) This value is used to identify a current J_1(V_1).
    
    (3) This current is used to evaluate an effective emission voltage V_em_1.
    
    (4) This emission voltage is then used to evaluate J_2(V - V_1).
    
    (5) This in turn is used to evaluate the effective emission voltage of the bottom cell (V_em_2)"""

    try:
        
        V[ 0 ]
        
    except:
        
        pass
    
    else:
        
        V = float( V[ 0 ] )
        
    #--------------------------------------------------------------------------
    # Iterate to Find Parameters
    #--------------------------------------------------------------------------    
    
    Delta_J = 100 * Convergence_Criterion
    
    Iteration = 0
        
    while Delta_J > Convergence_Criterion:
        
        #----------------------------------------------------------------------
        # Estimate V_1 by setting V_2 = V - V_1 and calculating J_1(V_1) = J_2(V_2)        
        #----------------------------------------------------------------------        
        
        if Iteration > Max_Iterations:
            
            break
        
        else:
            
            if Iteration == 0:      # Start by neglecting couplings
                
                V_em_1 = 0
                V_em_2 = 0
                                           
            V_1_prime = MJ_V_1_Calculator( V,
                               
                               EQEs_1, 
                               EQEs_2, 
                               
                               Incident_Photon_Fluxes,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_1,
                               EQE_EL_2,                   
                               
                               kT,
                               
                               R_sh_1,
                               R_sh_2,
                                               
                               R_s_1,
                               R_s_2,
                               Area,
                               
                               Iteration = Iteration,
                               
                               V_em_1 = V_em_1,
                               V_em_2 = V_em_2,
                               
                               Coupling_1_to_2 = Coupling_1_to_2,
                               Coupling_2_to_1 = Coupling_2_to_1 )
                                
            #----------------------------------------------------------------------
            # Calculate J_1(V_1), assuming V_em_2 = 0 (Negligible Backwards-emission)
            #----------------------------------------------------------------------        
                                        
            Incident_Photon_Fluxes_1 = Incident_Photon_Fluxes.copy()
            Incident_Photon_Fluxes_2 = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_1_prime, kT ) ) * Incident_Photon_Fluxes.copy()
            
            J_1 = MJ_J_i_Calculator( V_1_prime, EQEs_1, Incident_Photon_Fluxes_1,
                               
                            Black_Body_Fluxes, Energies, Wavelengths,
                               
                            EQE_EL_1, kT, R_sh_1,
                            
                            Coupling = Coupling_2_to_1,
                            EQEs_em = EQEs_2,
                            V_em = 0 )
            
            #----------------------------------------------------------------------
            # Recalculate = V_2_prime using current
            #----------------------------------------------------------------------        
            
            V_2_prime = V - V_1_prime - ( J_1 * Area * ( R_s_1 + R_s_2 ) ) / 1000
                        
            Incident_Photon_Fluxes_1 = Incident_Photon_Fluxes.copy()
            Incident_Photon_Fluxes_2 = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_1_prime, kT ) ) * Incident_Photon_Fluxes.copy()
            
            J_1 = MJ_J_i_Calculator( V_1_prime, EQEs_1, Incident_Photon_Fluxes_1,
                               
                            Black_Body_Fluxes, Energies, Wavelengths,
                               
                            EQE_EL_1, kT, R_sh_1,
                            
                            Coupling = Coupling_2_to_1,
                            EQEs_em = EQEs_2,
                            V_em = V_2_prime )
            
            #----------------------------------------------------------------------
            # Calculate V_em_1 Using this Current
            #----------------------------------------------------------------------        
            
            V_em_1 = MJ_V_em_i_Calculator( J_1, EQEs_1, Black_Body_Fluxes,
                                   
                                   Energies, Wavelengths,
                                   
                                   EQE_EL_1, kT, R_sh_1 )
            
            #----------------------------------------------------------------------
            # Calculate J_2 Using this emission voltage and V_2
            #----------------------------------------------------------------------        
                    
            J_2 = MJ_J_i_Calculator( V_2_prime, EQEs_2, Incident_Photon_Fluxes_2,
                               
                            Black_Body_Fluxes, Energies, Wavelengths,
                               
                            EQE_EL_2, kT, R_sh_2,
                            
                            Coupling = Coupling_1_to_2,
                            EQEs_em = EQEs_1,
                            V_em = V_em_1 )
            
            #----------------------------------------------------------------------
            # Calculate V_em_2 Using this Current
            #----------------------------------------------------------------------        
                                            
            V_em_2 = MJ_V_em_i_Calculator( J_2, EQEs_2, Black_Body_Fluxes,
                                   
                                   Energies, Wavelengths,
                                   
                                   EQE_EL_2, kT, R_sh_2 )                
        
#            print( Iteration, V_1_prime, V_2_prime, J_1 )
            
            #----------------------------------------------------------------------
            # Update Parameters
            #----------------------------------------------------------------------        
                                
            if Iteration == 0:
                            
                Delta_J = 100 * Convergence_Criterion
                
            else:
                
                Delta_J = abs( J_1 - J_Previous ) / abs( J_1 )
                    
            J_Previous = J_1
    
            Iteration += 1
            
    return J_1

def MJ_V_oc_i_Finder( EQEs_i, 
                   
                   Incident_Photon_Fluxes_i,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_i,
                   
                   kT,
                   
                   R_sh_i,
                   
                   Coupling = False,
                   
                   EQEs_em = 1 ):
    
    """Find the V_oc across a singular cell using the incident photon fluxes,
    
    EQEs, R_sh, etc. The emission voltage goes to nought at open-circuit 
    
    conditions."""
    
    def Function_to_Find_Root_of( V ):
        
        return MJ_J_i_Calculator( V,
                           
                           EQEs_i, 
                           
                           Incident_Photon_Fluxes_i,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL_i,
                           
                           kT,
                           
                           R_sh_i,
                           
                           Coupling = Coupling,
                           
                           EQEs_em = EQEs_em,
                           
                           V_em = 0 )
    
    Max_V = 1
    Bound_Found = False
    
    while not Bound_Found:
        
        try:
            
            brentq( Function_to_Find_Root_of, 0, Max_V )
            
        except:
            
            Max_V += 0.2
            
        else:
            
            Bound_Found = True
                
    return brentq( Function_to_Find_Root_of, 0, Max_V )

def MJ_V_oc_Finder(
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                    
                   R_s_1,
                   R_s_2,
                   Area,
                   
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Min_V = 1,
                   Max_V = 3,
                   
                   Convergence_Criterion = 1e-6,
                   Max_Iterations = 300 ):
    
    """Find the total V_oc using the total current."""
    
    def Function_to_Find_Root_of( V ):
        
        return MJ_J_tot_Calculator( V,
                           
                           EQEs_1, 
                           EQEs_2, 
                           
                           Incident_Photon_Fluxes,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL_1,
                           EQE_EL_2,                   
                           
                           kT,
                           
                           R_sh_1,
                           R_sh_2,
                           
                           R_s_1,
                           R_s_2,
                           Area,
                          
                           Coupling_1_to_2 = Coupling_1_to_2,
                           Coupling_2_to_1 = Coupling_2_to_1,
                           
                           Convergence_Criterion = Convergence_Criterion,
                           Max_Iterations = Max_Iterations )
    
    Lower_Bound_Found = False
    Upper_Bound_Found = False
    
    # Bound_Found = False
    # while not Bound_Found:
        
    #     try:
            
    #         brentq( Function_to_Find_Root_of, Min_V, Max_V )
            
    #     except:
            
    #         Max_V += 0.2
            
    #     else:
            
    #         Bound_Found = True
                
    # return brentq( Function_to_Find_Root_of, Min_V, Max_V )
    
    while not Lower_Bound_Found:
                    
        if Function_to_Find_Root_of( Min_V ) < 0:
            
            Lower_Bound_Found = True
            
        else:
            
            Min_V = 0.9 * Min_V 

    while not Upper_Bound_Found:
                    
        if Function_to_Find_Root_of( Max_V ) > 0:
            
            Upper_Bound_Found = True
            
        else:
            
            Max_V = 1.1 * Max_V
                        
    # USE BINARY SEARCH ALGORTIHM TO FIND Voc
            
    V_lower = Min_V
    V_upper = Max_V
    Delta_V = 1
    V_guess = ( V_lower + V_upper ) / 2
    
    if Lower_Bound_Found and Upper_Bound_Found:
        
        while Delta_V > Convergence_Criterion:
            
            if Function_to_Find_Root_of( V_guess ) > 0:
                
                V_upper = V_guess
                
            else: 
                
                V_lower = V_guess
                
            New_V_guess = ( V_lower + V_upper ) / 2

            Delta_V = abs( ( New_V_guess - V_guess ) / V_guess )
            
            V_guess = New_V_guess
            
        return V_guess
        
def MJ_V_mpp_Finder( 
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                    
                   R_s_1,
                   R_s_2,
                   Area,
                   
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Min_V = 1,
                   Max_V = 5,
                   
                   Convergence_Criterion = 1e-6,
                   Max_Iterations = 300 ):
    
    def Function_to_Minimize( V ):
        
        return V * MJ_J_tot_Calculator( V,
                           
                           EQEs_1, 
                           EQEs_2, 
                           
                           Incident_Photon_Fluxes,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL_1,
                           EQE_EL_2,                   
                           
                           kT,
                           
                           R_sh_1,
                           R_sh_2,
                           
                           R_s_1,
                           R_s_2,
                           Area,
                          
                           Coupling_1_to_2 = Coupling_1_to_2,
                           Coupling_2_to_1 = Coupling_2_to_1,
                           
                           Convergence_Criterion = Convergence_Criterion,
                           Max_Iterations = Max_Iterations )

    x0 = 0.8 * Max_V
            
    return minimize( Function_to_Minimize, x0, bounds = ( ( Min_V, Max_V ), ) )[ 'x' ][ 0 ]

def MJ_Two_Terminal_Figure_of_Merit_Calculator( Light_Power,
        
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                
                   R_s_1,
                   R_s_2,
                   Area,
                                      
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Convergence_Criterion = 1e-6,
                   Max_Iterations = 300 ):
    
    """Evaluate the photovoltaic figures-of-merit for a multi-junction device,
    
    with or without radiative couplings (from junction 1/top junction to 
                                         
    junction 2/bottom junction) using the EQEs of each cell, the incident 
    
    photon fluxes (in 10^-3 photons/cm2/s), the wavelengths (in), and the 
    
    thermal energy (kT, in eV)."""
    
    #--------------------------------------------------------------------------
    # Find Individual V_ocs
    #--------------------------------------------------------------------------    
        
    Incident_Photon_Fluxes_1_oc = Incident_Photon_Fluxes.copy()
    
    V_oc_1 = MJ_V_oc_i_Finder( EQEs_1, 
                       
                       Incident_Photon_Fluxes_1_oc,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_1,
                       
                       kT,
                       
                       R_sh_1,
                       
                       Coupling = Coupling_2_to_1,
                       
                       EQEs_em = EQEs_2  )
    
    Incident_Photon_Fluxes_2_oc = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_oc_1, kT ) ) * Incident_Photon_Fluxes.copy()
    
    V_oc_2 = MJ_V_oc_i_Finder( EQEs_2, 
                       
                       Incident_Photon_Fluxes_2_oc,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_2,
                       
                       kT,
                       
                       R_sh_2,
                       
                       Coupling = Coupling_1_to_2,
                       
                       EQEs_em = EQEs_1 )
    
    #--------------------------------------------------------------------------
    # Find Total V_oc
    #--------------------------------------------------------------------------    
    
    # Define voltage region-of-interest:
        
    Min_V = max( V_oc_1, V_oc_2 )
    V_oc_sum = V_oc_1 + V_oc_2    
    
#    print( 'Vocs: ', V_oc_1, V_oc_2, Min_V, V_oc_sum )
    
    V_oc = MJ_V_oc_Finder(
                    EQEs_1, 
                    EQEs_2, 
               
                    Incident_Photon_Fluxes,
               
                    Black_Body_Fluxes,
               
                    Energies,
               
                    Wavelengths,
               
                    EQE_EL_1,
                    EQE_EL_2,                   
               
                    kT,
               
                    R_sh_1,
                    R_sh_2,
                
                    R_s_1,
                    R_s_2,
                    Area,
               
                    Coupling_1_to_2 = Coupling_1_to_2,
                    Coupling_2_to_1 = Coupling_2_to_1,
               
                    Min_V = 1.2 * Min_V,
                    Max_V = V_oc_sum, 
               
                    Convergence_Criterion = Convergence_Criterion,
                    Max_Iterations = Max_Iterations )
    
    #--------------------------------------------------------------------------
    # Find V_mpp
    #--------------------------------------------------------------------------    
    
    Max_V = V_oc
    
    V_mpp = MJ_V_mpp_Finder( 
                       EQEs_1, 
                       EQEs_2, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_1,
                       EQE_EL_2,                   
                       
                       kT,
                       
                       R_sh_1,
                       R_sh_2,
                        
                       R_s_1,
                       R_s_2,
                       Area,
                       
                       Coupling_1_to_2 = Coupling_1_to_2,
                       Coupling_2_to_1 = Coupling_2_to_1,
                       
                       Min_V = Min_V,
                       Max_V = V_oc_sum,
                       
                       Convergence_Criterion = Convergence_Criterion,
                       Max_Iterations = Max_Iterations )
            
    #--------------------------------------------------------------------------
    # Find J_mpp
    #--------------------------------------------------------------------------    
    
    J_mpp = MJ_J_tot_Calculator( V_mpp,
                       
                       EQEs_1, 
                       EQEs_2, 
                       
                       Incident_Photon_Fluxes,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_1,
                       EQE_EL_2,                   
                       
                       kT,
                       
                       R_sh_1,
                       R_sh_2,
                       
                       R_s_1,
                       R_s_2,
                       Area,
                      
                       Coupling_1_to_2 = Coupling_1_to_2,
                       Coupling_2_to_1 = Coupling_2_to_1,
                       
                       Convergence_Criterion = Convergence_Criterion,
                       Max_Iterations = Max_Iterations )
    
    #--------------------------------------------------------------------------
    # Calculate and Return Remaining Figures-of-Merit
    #--------------------------------------------------------------------------    
        
    P_mpp = abs( V_mpp * J_mpp ) # mW/cm2
    
    PCE = P_mpp / Light_Power
        
    return { "V_oc" : V_oc,
            
            "V_mpp" : V_mpp,
            
             "J_mpp" : J_mpp,
             
             "P_mpp" : P_mpp,
             
             "PCE" : PCE }

def MJ_Lone_Cell_V_mpp_Calculator( EQEs_i, 
    
    Incident_Photon_Fluxes_i,
    
    Black_Body_Fluxes,
    
    Energies,
    
    Wavelengths,
    
    EQE_EL_i,
    
    kT,
    
    R_sh_i,
    
    max_V,
    
    min_V = 0,
 
    Coupling = False,
    
    EQEs_em = 1,
    
    V_em = 0 ):
 
    """Evaluate the maximum power point voltage of a lone cell using
    
    its EQEs and incident fluxes."""

    def Function_to_Minimise( V_i ):
        
        return V_i * MJ_J_i_Calculator( V_i,
                           
                           EQEs_i, 
                           
                           Incident_Photon_Fluxes_i,
                           
                           Black_Body_Fluxes,
                           
                           Energies,
                           
                           Wavelengths,
                           
                           EQE_EL_i,
                           
                           kT,
                           
                           R_sh_i,
                           
                           Coupling = Coupling,
                           
                           EQEs_em = EQEs_em,
                           
                           V_em = V_em )
        
    x0 = 0.8 * max_V
    
    return minimize( Function_to_Minimise, x0, bounds = ( ( min_V, max_V ), ) )[ 'x' ][ 0 ]
    
def MJ_Four_Terminal_V_mpps_Calculator( 
        
        EQEs_1, 
        EQEs_2, 
        
        Incident_Photon_Fluxes,
        
        Black_Body_Fluxes,
        
        Energies,
        
        Wavelengths,
        
        EQE_EL_1,
        EQE_EL_2,                   
        
        kT,
        
        R_sh_1,
        R_sh_2,
        
        V_oc_1,
        V_oc_2,
                
        Coupling_1_to_2 = False,
        Coupling_2_to_1 = False,
        
        Convergence_Criterion = 1e-6,
        Max_Iterations = 300 ):
    
    """Determine the maximum power point voltages of a two-junction device in
    
    the four-terminal configuration using an iterative process.
    
    Assume in the first iteration that the maximum power point parameters can 
    
    be determined without radiative coupling, then introduce these effects
    
    as a second-order process."""
    
    #--------------------------------------------------------------------------
    # First, Calculate V_mpp of Cells in Isolation (No Coupling)
    #--------------------------------------------------------------------------

    Incident_Photon_Fluxes_1_mpp = Incident_Photon_Fluxes.copy()
    
    V_mpp_1_Lone = MJ_Lone_Cell_V_mpp_Calculator( EQEs_1, 
        
        Incident_Photon_Fluxes_1_mpp,
        
        Black_Body_Fluxes,
        
        Energies,
        
        Wavelengths,
        
        EQE_EL_1,
        
        kT,
        
        R_sh_1,
        
        V_oc_1,
        
        Coupling = False,   ##### ASSUME COUPLING IS FALSE TO BEGIN
        
        EQEs_em = 1,
        
        V_em = 0 )
    
    Incident_Photon_Fluxes_2_mpp = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_mpp_1_Lone, kT ) ) * Incident_Photon_Fluxes.copy()

    V_mpp_2_Lone = MJ_Lone_Cell_V_mpp_Calculator( EQEs_2, 
        
        Incident_Photon_Fluxes_2_mpp,
        
        Black_Body_Fluxes,
        
        Energies,
        
        Wavelengths,
        
        EQE_EL_2,
        
        kT,
        
        R_sh_2,
        
        V_oc_2,
        
        Coupling = False,   ##### ASSUME COUPLING IS FALSE TO BEGIN
        
        EQEs_em = 1,
        
        V_em = 0 ) 
    
    J_mpp_1_Lone = MJ_J_i_Calculator( V_mpp_1_Lone,
           
           EQEs_1, 
           
           Incident_Photon_Fluxes_1_mpp,
           
           Black_Body_Fluxes,
           
           Energies,
           
           Wavelengths,
           
           EQE_EL_1,
           
           kT,
           
           R_sh_1,
           
           Coupling = False,
           
           EQEs_em = 1,
           
           V_em = 0)

    J_mpp_2_Lone = MJ_J_i_Calculator( V_mpp_2_Lone,
           
           EQEs_2, 
           
           Incident_Photon_Fluxes_2_mpp,
           
           Black_Body_Fluxes,
           
           Energies,
           
           Wavelengths,
           
           EQE_EL_2,
           
           kT,
           
           R_sh_2,
           
           Coupling = False,
           
           EQEs_em = 1,
           
           V_em = 0)
    
    #--------------------------------------------------------------------------
    # If No Coupling, Return Result
    #--------------------------------------------------------------------------
    
    if Coupling_1_to_2 == False and Coupling_2_to_1 == False:
        
        return V_mpp_1_Lone, V_mpp_2_Lone, J_mpp_1_Lone, J_mpp_2_Lone
    
    else:
        
        #----------------------------------------------------------------------
        # Iterate Until Converged
        #----------------------------------------------------------------------
        
        Iteration = 0 
        Delta_P = 100 * Convergence_Criterion
        
        while Delta_P > Convergence_Criterion:
                    
            #------------------------------------------------------------------
            # Evaluate J_mpp of both Cells
            #------------------------------------------------------------------
            
            if Iteration == 0:
                
                #--------------------------------------------------------------
                # Assume No Coupling in Zeroth Iteration
                #--------------------------------------------------------------
                
                J_mpp_1 = MJ_J_i_Calculator( V_mpp_1_Lone,
                                   
                                   EQEs_1, 
                                   
                                   Incident_Photon_Fluxes_1_mpp,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_1,
                                   
                                   kT,
                                   
                                   R_sh_1,
                                   
                                   Coupling = False,
                                   
                                   EQEs_em = 1,
                                   
                                   V_em = 0 )
                
                J_mpp_2 = MJ_J_i_Calculator( V_mpp_2_Lone,
                                   
                                   EQEs_2, 
                                   
                                   Incident_Photon_Fluxes_2_mpp,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_2,
                                   
                                   kT,
                                   
                                   R_sh_2,
                                   
                                   Coupling = False,
                                   
                                   EQEs_em = 1,
                                   
                                   V_em = 0 )
                
            else:  
                
                #--------------------------------------------------------------
                # Include Coupling in Other Iterations
                #--------------------------------------------------------------
                                    
                J_mpp_1 = MJ_J_i_Calculator( V_mpp_1,
                                   
                                   EQEs_1, 
                                   
                                   Incident_Photon_Fluxes_1_mpp,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_1,
                                   
                                   kT,
                                   
                                   R_sh_1,
                                   
                                   Coupling = Coupling_2_to_1,
                                   
                                   EQEs_em = EQEs_2,
                                   
                                   V_em = V_em_2 )
                
                J_mpp_2 = MJ_J_i_Calculator( V_mpp_2,
                                   
                                   EQEs_2, 
                                   
                                   Incident_Photon_Fluxes_2_mpp,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_2,
                                   
                                   kT,
                                   
                                   R_sh_2,
                                   
                                   Coupling = Coupling_1_to_2,
                                   
                                   EQEs_em = EQEs_1,
                                   
                                   V_em = V_em_1 )
            
            #------------------------------------------------------------------
            # Evaluate Emission Voltages of Both Cells Using Reciprocity
            #------------------------------------------------------------------
                            
            V_em_1 = MJ_V_em_i_Calculator( J_mpp_1,
                                                       
                                       EQEs_1, 
                                                              
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
                                       
                                       EQE_EL_1,
                                       
                                       kT,
                                       
                                       R_sh_1 )
            
            V_em_2 = MJ_V_em_i_Calculator( J_mpp_2,
                                                       
                                       EQEs_2, 
                                                              
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
                                       
                                       EQE_EL_2,
                                       
                                       kT,
                                       
                                       R_sh_2 )
            
            #------------------------------------------------------------------
            # Re-Evaluate Maximum Power Point Voltages
            #------------------------------------------------------------------
    
            Incident_Photon_Fluxes_1_mpp = Incident_Photon_Fluxes.copy()
            
            V_mpp_1 = MJ_Lone_Cell_V_mpp_Calculator( EQEs_1, 
                
                Incident_Photon_Fluxes_1_mpp,
                
                Black_Body_Fluxes,
                
                Energies,
                
                Wavelengths,
                
                EQE_EL_1,
                
                kT,
                
                R_sh_1,
                
                V_oc_1,
                
                Coupling = Coupling_2_to_1,
                
                EQEs_em = EQEs_2,
                
                V_em = V_em_2 )
            
            Incident_Photon_Fluxes_2_mpp = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_mpp_1, kT ) ) * Incident_Photon_Fluxes.copy()
        
            V_mpp_2 = MJ_Lone_Cell_V_mpp_Calculator( EQEs_2, 
                
                Incident_Photon_Fluxes_2_mpp,
                
                Black_Body_Fluxes,
                
                Energies,
                
                Wavelengths,
                
                EQE_EL_2,
                
                kT,
                
                R_sh_2,
                
                V_oc_2,
                
                Coupling = Coupling_1_to_2,   ##### ASSUME COUPLING IS FALSE TO BEGIN
                
                EQEs_em = EQEs_1,
                
                V_em = V_em_1 )
            
            #------------------------------------------------------------------
            # Evaluate Maximum Power Point Power
            #------------------------------------------------------------------
    
            J_mpp_1 = MJ_J_i_Calculator( V_mpp_1,
                   
                   EQEs_1, 
                   
                   Incident_Photon_Fluxes_1_mpp,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   
                   kT,
                   
                   R_sh_1,
                   
                   Coupling = Coupling_2_to_1,
                   
                   EQEs_em = EQEs_2,
                   
                   V_em = V_em_2 )

            J_mpp_2 = MJ_J_i_Calculator( V_mpp_2,
                   
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes_2_mpp,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_2,
                   
                   kT,
                   
                   R_sh_2,
                   
                   Coupling = Coupling_1_to_2,
                   
                   EQEs_em = EQEs_1,
                   
                   V_em = V_em_1 )
        
            P_mpp_1 = abs( J_mpp_1 ) * V_mpp_1
            P_mpp_2 = abs( J_mpp_2 ) * V_mpp_2
            
            if Iteration == 0:
                    
                Delta_P = 100 * Convergence_Criterion
                
                P_mpp_tot_Prior = P_mpp_1 + P_mpp_2
            
            else:
                
                P_mpp_tot = P_mpp_1 + P_mpp_2 
                
                Delta_P = abs( P_mpp_tot - P_mpp_tot_Prior ) / P_mpp_tot_Prior
                
                P_mpp_tot_Prior = P_mpp_tot
                    
            # print( Iteration, P_mpp_1, P_mpp_2, P_mpp_tot_Prior, Delta_P )
            
            Iteration += 1
    
        return V_mpp_1, V_mpp_2, J_mpp_1, J_mpp_2
            
def MJ_Four_Terminal_Figure_of_Merit_Calculator( Light_Power,
        
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                                      
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Convergence_Criterion = 1e-6,
                   Max_Iterations = 300 ):
    
    """Evaluate the photovoltaic figures-of-merit for a multi-junction device
    
    in the four-terminal case (where the two junctions have seperate outputs),
    
    with or without radiative couplings (from junction 1/top junction to 
                                         
    junction 2/bottom junction) using the EQEs of each cell, the incident 
    
    photon fluxes (in 10^-3 photons/cm2/s), the wavelengths (in), and the 
    
    thermal energy (kT, in eV).
    
    
    The function is broken into two parts: in the first , band-filling is False,
    
    and the computation time is reduced. In the second case, band-filling
    
    is True."""        
    
    #--------------------------------------------------------------------------
    # Find Individual V_ocs
    #--------------------------------------------------------------------------    
        
    Incident_Photon_Fluxes_1_oc = Incident_Photon_Fluxes.copy()
    
    V_oc_1 = MJ_V_oc_i_Finder( EQEs_1, 
                       
                       Incident_Photon_Fluxes_1_oc,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_1,
                       
                       kT,
                       
                       R_sh_1,
                       
                       Coupling = Coupling_2_to_1,
                       
                       EQEs_em = EQEs_2  )
    
    Incident_Photon_Fluxes_2_oc = ( 1 - EQEs_1 * Band_Filling_Factor( Energies, V_oc_1, kT ) ) * Incident_Photon_Fluxes.copy()
    
    V_oc_2 = MJ_V_oc_i_Finder( EQEs_2, 
                       
                       Incident_Photon_Fluxes_2_oc,
                       
                       Black_Body_Fluxes,
                       
                       Energies,
                       
                       Wavelengths,
                       
                       EQE_EL_2,
                       
                       kT,
                       
                       R_sh_2,
                       
                       Coupling = Coupling_1_to_2,
                       
                       EQEs_em = EQEs_1 )
    
    #--------------------------------------------------------------------------
    # Find V_mpp of Each Cell in Isolation
    #--------------------------------------------------------------------------    
    
    V_mpp_1, V_mpp_2, J_mpp_1, J_mpp_2 = MJ_Four_Terminal_V_mpps_Calculator( 
            
            EQEs_1, 
            EQEs_2, 
            
            Incident_Photon_Fluxes,
            
            Black_Body_Fluxes,
            
            Energies,
            
            Wavelengths,
            
            EQE_EL_1,
            EQE_EL_2,                   
            
            kT,
            
            R_sh_1,
            R_sh_2,
            
            V_oc_1,
            V_oc_2,
                    
            Coupling_1_to_2 = Coupling_1_to_2,
            Coupling_2_to_1 = Coupling_2_to_1,
            
            Convergence_Criterion = Convergence_Criterion,
            Max_Iterations = Max_Iterations )
    
    #--------------------------------------------------------------------------
    # Find P_mpps and PCE
    #--------------------------------------------------------------------------    
    
    P_mpp_1 = abs( J_mpp_1 ) * V_mpp_1
    P_mpp_2 = abs( J_mpp_2 ) * V_mpp_2
    
    P_mpp = P_mpp_1 + P_mpp_2

    PCE = P_mpp / Light_Power
        
    return { "V_oc_1" : V_oc_1,
            
            "V_oc_2" : V_oc_2,
            
            "V_mpp_1" : V_mpp_1,
            
            "V_mpp_2" : V_mpp_2,
            
             "J_mpp_1" : J_mpp_1,

             "J_mpp_2" : J_mpp_2,
             
             "P_mpp_1" : P_mpp_1,

             "P_mpp_2" : P_mpp_2,
             
             "P_mpp" : P_mpp,

             "PCE" : PCE }

def MJ_Small_Area_Figure_of_Merit_Calculator( Light_Power,
        
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,

                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   
                   kT,
                   
                   
                   R_sh_1,
                   R_sh_2,
                                
                   R_s_1,
                   R_s_2,
                   Area_celL,

                   Terminal_Type = "Two-Terminal (Serially-Connected)",
                   
                   Coupling_1_to_2 = False,
                   Coupling_2_to_1 = False,
                   
                   Convergence_Criterion = 1e-6, 
                   Max_Iterations = 300 ):
    

    if Terminal_Type == "Two-Terminal (Serially-Connected)":
        
        return MJ_Two_Terminal_Figure_of_Merit_Calculator( Light_Power,
         
                    EQEs_1, 
                    EQEs_2, 
                    
                    Incident_Photon_Fluxes,
                    
                    Black_Body_Fluxes,
                    
                    Energies,
                    
                    Wavelengths,
                    
                    EQE_EL_1,
                    EQE_EL_2,                   
                    
                    kT,
                    
                    R_sh_1,
                    R_sh_2,
                 
                    R_s_1,
                    R_s_2,
                    Area_celL,
                                      
                    Coupling_1_to_2 = Coupling_1_to_2,
                    Coupling_2_to_1 = Coupling_2_to_1,
                   
                    Convergence_Criterion = Convergence_Criterion,
                    Max_Iterations = Max_Iterations )

    elif Terminal_Type == "Four-Terminal":

        return MJ_Four_Terminal_Figure_of_Merit_Calculator( Light_Power,
        
                   EQEs_1, 
                   EQEs_2, 
                   
                   Incident_Photon_Fluxes,
                   
                   Black_Body_Fluxes,
                   
                   Energies,
                   
                   Wavelengths,
                   
                   EQE_EL_1,
                   EQE_EL_2,                   
                   
                   kT,
                   
                   R_sh_1,
                   R_sh_2,
                                      
                   Coupling_1_to_2 = Coupling_1_to_2,
                   Coupling_2_to_1 = Coupling_2_to_1,
                   
                   Convergence_Criterion = Convergence_Criterion,
                   Max_Iterations = Max_Iterations )

#%%===========================================================================7
# 7. Scaling-Up Devices and Systems
#==============================================================================

def Large_Area_FoM_Calculator( Small_Area_Output, I_light,
                              
                   Area_Module, Area_cell, N_cells,
                   
                   L_x, L_y, delta,
                                      
                   Evaluate_System_Efficiency = True ):
    
    """I_light is irradiance in mW/cm2.
        
    delta is interconnect width in cm.
    
    Area of the module of N cells is in cm2.
    
    R_sh in Ohm cm2,
    
    R_sheet in Ohm / m2,
    
    R_s in Ohm,
    
    Whether or not the system efficiency of an energy-harvesting board
    
    is to be evaluated is given as True or False."""
    
    #--------------------------------------------------------------------------    
    # Determine Incident Light Power and Area of Each Cell    
    #--------------------------------------------------------------------------
    
    P_light = I_light * Area_Module # units are mW
        
    #--------------------------------------------------------------------------    
    # Determine Cell Output
    #--------------------------------------------------------------------------
    
    if Area_cell > 0:
        
        Output = Small_Area_Output
        
    else:
        
        Output = { "J_sc" : 0,
                
                "V_oc" : 0,
                
                "V_mpp" : 0,
                
                 "J_mpp" : 0,
                 
                 "P_mpp" : 0,
                 
                 "PCE" : 0 } 
    
    #--------------------------------------------------------------------------    
    # Determine Output Voltages
    #--------------------------------------------------------------------------
        
    V_oc_Module = Output[ 'V_oc' ] * N_cells
    
    V_mpp_Module = Output[ 'V_mpp' ] * N_cells

    #--------------------------------------------------------------------------    
    # Determine Output Current (Not Current Density)
    #--------------------------------------------------------------------------
                
    J_sc = abs( Output[ 'J_sc' ] ) # Units are mA/cm2
    
    J_mpp = Output[ 'J_mpp' ] # Units are mA/cm2
    
    I_sc_Module = J_sc * Area_cell # Units are mA; Assume current same in all cells
    
    I_mpp_Module = J_mpp * Area_cell # Units are mA; Assume current same in all cells
    
    #--------------------------------------------------------------------------    
    # Determine Remaining FoMs
    #--------------------------------------------------------------------------
          
    P_mpp_Module = abs( I_mpp_Module * V_mpp_Module ) # Units are mW
    
    PCE_Module = P_mpp_Module / P_light
    
    Output = Output | { 'Area_Module' : Area_Module,
                       
            'N_cells' : N_cells,
            
            'Cell_Width' : ( L_x - ( N_cells - 1 ) * delta ) / N_cells,

            'Area_cell' : Area_cell,
                        
            'Area_interconnects' : ( N_cells - 1 ) * delta * L_y,
            
            'I_sc_Module' : I_sc_Module,
                        
             'V_oc_Module' : V_oc_Module,

             'V_mpp_Module' : V_mpp_Module,
 
             'I_mpp_Module' : I_mpp_Module,
             
             'P_mpp_Module' : P_mpp_Module,
             
             'PCE_Module' : PCE_Module }
    
    if Evaluate_System_Efficiency:
                
        eta_100uA, eta_1mA, eta_10mA, eta_interp = Board_Efficiency_Calculator( V_mpp_Module, abs( I_mpp_Module ) )

        Output = Output | { 'eta_interp' : eta_interp,
                           
                           'eta_system' : eta_interp * PCE_Module,
                           
                           'eta_100uA' : eta_100uA,
                           
                           'eta_1mA' : eta_1mA,
                           
                           'eta_10mA' : eta_10mA }
        
    return Output

#%%===========================================================================8
# 8. Generating a GUI
#==============================================================================

"""
In this section, the widgets for the GUI are compiled. The descriptions, default values,

current values, and 

"""

# Write Path to Widget Inputs File

Parameter_Filepath = path.join( Current_Working_Directory, "Parameter_Values.xlsx" )

Widgets = {}

#------------------------------------------------------------------------------
# Determine Whether GUI Must Be Initialised (i.e., is "GUI_Initialiser.py" run)
#------------------------------------------------------------------------------

GUI = False

if "Initialise_GUI.txt" in listdir( Current_Working_Directory ):
    
    GUI = True
    
    remove( path.join( Current_Working_Directory, "Initialise_GUI.txt" ) )
    
if GUI:
        
    #------------------------------------------------------------------------------
    # Create Dropdown Widgets
    #------------------------------------------------------------------------------
    
    Dropdown_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Dropdowns', engine='openpyxl' )
    
    for j in range( len( Dropdown_Widget_DF ) ):
        
        Row = Dropdown_Widget_DF.iloc[ j ]
        
        Description = Row[ "Description" ] 
        
        #==========================================================================
        # Add light sources using available data
        #==========================================================================
        
        if Description == 'Light Source':
            
            if str( Row[ "Saved Value" ] ) in Light_Sources:
                
                Value = str( Row[ "Saved Value" ] )
                
            else:
                
                Value = 'AM1.5G'
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = Value,
                                                       
                                                       options = Light_Sources,
                                                   
                                                       style = { 'description_width' : 'initial' } )
            
        #==========================================================================
        # Add EQE options using available data
        #==========================================================================
        
        elif Description == 'Experimental EQE Model':
            
            if str( Row[ "Saved Value" ] ) in System_Identifiers:
                
                Value = str( Row[ "Saved Value" ] )
                
            else:
                
                Value = 'INO_1'
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = Value,
                                                       
                                                       options = System_Identifiers,
                                                   
                                                       style = { 'description_width' : 'initial' } )        
    
        #==========================================================================
        # Add Bottom Cell EQE options using available data
        #==========================================================================
        
        elif Description == 'Bottom Cell Experimental EQE Model':
            
            if str( Row[ "Saved Value" ] ) in System_Identifiers:
                
                Value = str( Row[ "Saved Value" ] )
                
            else:
                
                Value = 'INO_1'
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = Value,
                                                       
                                                       options = System_Identifiers,
                                                   
                                                       style = { 'description_width' : 'initial' } )   
            
        #==========================================================================
        # Add Material Options for Optical Modelling Widgets
        #==========================================================================
            
        elif Description == 'Surrounding Medium':
            
            if str( Row[ "Saved Value" ] ) in Materials:
                
                Value = str( Row[ "Saved Value" ] )
                
            else:
                
                Value = 'air'
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = Value,
                                                       
                                                       options = Materials,
                                                   
                                                       style = { 'description_width' : 'initial' } )   
            
        elif 'Layer Type ' in Description:
                    
            if str( Row[ "Saved Value" ] ) in Materials:
                
                Value = str( Row[ "Saved Value" ] )
                
            else:
                
                Value = 'glass'
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = Value,
                                                       
                                                       options = Materials,
                                                   
                                                       style = { 'description_width' : 'initial' } )   
                    
        #==========================================================================
        # Otherwise, Use Options in Excel Sheet
        #==========================================================================
                
        else: 
            
            Widgets[ Description ] = Dropdown( description = Description + ': ',
                                                       
                                                       value = str( Row[ "Saved Value" ] ),
                                                       
                                                       options = Row[ "Options" ].split( ', ' ),
                                                   
                                                       style = { 'description_width' : 'initial' } )
    
    #------------------------------------------------------------------------------
    # Create Float Widgets
    #------------------------------------------------------------------------------

    Float_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Floats', engine='openpyxl' )
    
    for j in range( len( Float_Widget_DF ) ):
        
        Row = Float_Widget_DF.iloc[ j ]
        
        Widgets[ Row[ "Parameter" ] ] = BoundedFloatText( description = str( Row[ "Description" ] ) + ': ',
                                                   
                                                   value = Row[ "Saved Value" ],
                                                   
                                                   min = Row[ "Minimum Value" ],
                                                   
                                                   max = Row[ "Maximum Value" ],
                                               
                                                   style = { 'description_width' : 'initial' } )
    
    #------------------------------------------------------------------------------
    # Create Int Widgets
    #------------------------------------------------------------------------------

    Int_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Integers', engine='openpyxl')
    
    for j in range( len( Int_Widget_DF ) ):
        
        Row = Int_Widget_DF.iloc[ j ]
        
        Widgets[ Row[ "Parameter" ] ] = BoundedIntText( description = str( Row[ "Description" ] ) + ': ',
                                                   
                                                   value = Row[ "Saved Value" ],
                                                   
                                                   min = Row[ "Minimum Value" ],
                                                   
                                                   max = Row[ "Maximum Value" ],
                                               
                                                   style = { 'description_width' : 'initial' } )    
        
    
    #==========================================================================
    # GUI Control Function
    #==========================================================================
    
    def GUI_Updater( Change ):
    
        GUI_Mode = Widgets[ "Simulation Mode" ].value
        
        #----------------------------------------------------------------------
        # Hide All Widgets But Unchanging Ones
        #----------------------------------------------------------------------
        
        Dont_Hide_Widgets = ["Simulation Mode",
                             'Optically Modelled',
                             "Junction Type", 
                             "Module Simulation",
                             'Generate Plots',
                             'R_series (Ohms)',
                             'R_shunt (Ohms cm2)',
                             "Include Point Defects",
                             'Intensity Unit',
                             'Minimum Energy (eV)',
                             'Maximum Energy (eV)',
                             'Number of Energies',
                             'Temperature',
                             'Temperature Unit'
                             ]
    
        for Widget_Name, Widget in Widgets.items():
            
            if Widget_Name not in Dont_Hide_Widgets:
                
                Widget.layout.display = 'none'
                
            else:
                
                Widget.layout.display = None
    
        Top_Cell_Label.layout.display = 'none'
        Bottom_Cell_Label.layout.display = 'none'
        Gap_Inputs_Label.layout.display = 'none'
                
        #----------------------------------------------------------------------
        # Reveal Module Widgets
        #----------------------------------------------------------------------
                
        if Widgets[ 'Module Simulation' ].value == 'True':
            
            Widgets[ 'Number of Cells' ].layout.display = None
            Widgets[ 'Module Width (cm)' ].layout.display = None
            Widgets[ 'Module Length (cm)'].layout.display = None
            Widgets[ 'Interconnect Width (cm)' ].layout.display = None
            
        else:
            
            Widgets[ 'Area (cm2)' ].layout.display = None
    
        #----------------------------------------------------------------------
        # If Mode is Has Varying Material, Do Not Allow Optical Modelling
        #----------------------------------------------------------------------
        
        if GUI_Mode not in [ 'Single Simulation', 'Varied Intensity', 
                        
                        'Varied Light Source', 'Typical Indoor Day' ]:
            
            Widgets[ 'Optically Modelled' ].value = 'False'
            Widgets[ 'Optically Modelled' ].layout.display = 'none'
    
        #----------------------------------------------------------------------
        # Reveal Material Widgets
        #----------------------------------------------------------------------
        
        Widgets[ 'Non-Radiative Loss Mode' ].layout.display = None
        
        if Widgets[ 'Optically Modelled' ].value != 'True':
            
            #----------------------------------------------------------------------
            # If not optically modelled, simulate or load an EQE:
            #----------------------------------------------------------------------            
                
            Widgets[ "Material Model" ].layout.display = None
            
            if Widgets[ "Material Model" ].value == 'Simulated':
                
                Simualted_EQE_Mode = Widgets[ "Simulated EQE Model" ].value
                Widgets[ "Simulated EQE Model" ].layout.display = None    
                
                if GUI_Mode != 'Varied Gap':
                    Widgets[ 'E_g (eV)' ].layout.display = None
                else:
                    Widgets[ 'Minimum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Maximum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Number of E_gs' ].layout.display = None     
                    
                Widgets[ 'Above Gap EQE' ].layout.display = None
                
                if Simualted_EQE_Mode == 'Step Function':
                    Widgets[ 'Below Gap EQE' ].layout.display = None
        
                elif Simualted_EQE_Mode == 'Urbach Tail':
                    Widgets[ 'Urbach Energy (eV)' ].layout.display = None
                    
                elif Simualted_EQE_Mode == 'Disordered Organic PV':
                    Widgets[ 'Excitonic Disorder (eV)' ].layout.display = None
                        
            else:
                
                Widgets[ "Experimental EQE Model" ].layout.display = None
    
        
        if Widgets[ 'Optically Modelled' ].value == 'True':
            
            #----------------------------------------------------------------------
            # If optically modelled, reveal inputs for simulating an EQE:
            #----------------------------------------------------------------------   
             
            Widgets[ 'Position Spacing (nm)' ].layout.display = None
            Widgets[ "Surrounding Medium" ].layout.display = None
            Widgets[ "Incident Angle (deg)" ].layout.display = None
            Widgets[ "Number of Layers in Device" ].layout.display = None
            Widgets[ "Incoherent Substrate" ].layout.display = None
            
            if Widgets[ "Incoherent Substrate" ].value == 'True':
                Widgets[ "Substrate Index" ].layout.display = None
    
            Widgets[ "Infinite Final Layer" ].layout.display = None
            
            Widgets[ "Active Layer Index" ].layout.display = None
            Widgets[ "Internal Quantum Efficiency" ].layout.display = None
    
            if Widgets[ "Junction Type" ].value == 'Multi-Junction':
                Widgets[ "Bottom Active Layer Index" ].layout.display = None
                Widgets[ "Bottom Cell Non-Radiative Loss Mode" ].layout.display = None
                
            for j in range( Widgets[ "Number of Layers in Device" ].value ):
                
                Widgets[ "Layer Type " + str( j + 1 ) ].layout.display = None
                Widgets[ "Layer Thickness " + str( j + 1 ) + " (nm)" ].layout.display = None
    
        #----------------------------------------------------------------------
        # Reveal Device Parameter Widgets
        #----------------------------------------------------------------------
        
        if Widgets[ "Include Point Defects" ].value == 'True':
            
            Widgets[ 'Point Defect Resistance (Ohms cm2)' ].layout.display = None
            Widgets[ "Simulate Number of Defects" ].layout.display = None
            if Widgets[ "Simulate Number of Defects" ].value == 'True':
                
                Widgets[ 'Point Defect Proclivity f_x (1/cm)' ].layout.display = None
                Widgets[ 'Point Defect Proclivity f_y (1/cm)' ].layout.display = None
                Widgets[ 'A_onset (cm)' ].layout.display = None
            else:
                Widgets[ "Number of Defects" ].layout.display = None
    
        #----------------------------------------------------------------------
        # Reveal Multi-Junction Widgets if True
        #----------------------------------------------------------------------
                
        if Widgets[ 'Junction Type' ].value == 'Multi-Junction':
            
            # Dont simulate/import EQE if Optically Modelled
            if Widgets[ 'Optically Modelled' ].value != 'True':
    
                Top_Cell_Label.layout.display = None
                Bottom_Cell_Label.layout.display = None
                
                Widgets[ "Multi-Junction Terminal Type" ].layout.display = None
                Widgets[ "Bottom Cell Material Model" ].layout.display = None
                Widgets[ "Bottom Cell Non-Radiative Loss Mode" ].layout.display = None
                Widgets[ 'Radiative Coupling Bottom to Top' ].layout.display = None
                Widgets[ 'Radiative Coupling Top to Bottom' ].layout.display = None
                Widgets[ 'Bottom Cell R_series (Ohms)' ].layout.display = None
                Widgets[ 'Bottom Cell R_shunt (Ohms cm2)' ].layout.display = None
                
                # Reveal Bottom Cell Input Widgets
                
                if Widgets[ "Bottom Cell Material Model" ].value == 'Simulated':
                    
        
                    Widgets[ "Bottom Cell Simulated EQE Model" ].layout.display = None
                    Widgets[ 'Bottom Cell E_g (eV)' ].layout.display = None
                    Widgets[ 'Bottom Cell Above Gap EQE' ].layout.display = None
                    
                    if Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Step Function':
                        
                        Widgets[ 'Bottom Cell Below Gap EQE' ].layout.display = None
        
                    elif Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Urbach Tail':
                        
                        Widgets[ 'Bottom Cell Urbach Energy (eV)' ].layout.display = None
                        
                    elif Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Disordered Organic PV':
                        
                        Widgets[ 'Bottom Cell Excitonic Disorder (eV)' ].layout.display = None
                            
                else:
                    
                    Widgets[ "Bottom Cell Experimental EQE Model" ].layout.display = None
    
        #----------------------------------------------------------------------
        # Reveal Remaining Widgets as Necessary 
        #----------------------------------------------------------------------
        
        if GUI_Mode in [ 'Single Simulation', 'Many-System Comparison', 'Varied Gap', 
                        
                        'Varied Lower Limit', "Multi-Junction - One Fixed Cell" ]:
            
            Widgets[ 'Light Source' ].layout.display = None
            Widgets[ 'Light Intensity' ].layout.display = None
            
        if GUI_Mode == 'Varied Gap':
            Widgets[ 'Material Model' ].value = 'Simulated'
            Widgets[ 'Material Model' ].layout.display = 'none'
            Widgets[ "Junction Type" ].value = 'Single-Junction' # Don't allow MJ in this case
            Widgets[ "Junction Type" ].layout.display = 'none'   
            
        if GUI_Mode == 'Varied Lower Limit':
            Widgets[ 'Material Model' ].value = 'Experimental (Loaded)'
            Widgets[ 'Material Model' ].layout.display = 'none'
            Widgets[ "Junction Type" ].value = 'Single-Junction' # Don't allow MJ in this case
            Widgets[ "Junction Type" ].layout.display = 'none'           
            
        if GUI_Mode == 'Varied Intensity':
            
            Widgets[ 'Light Source' ].layout.display = None
            Widgets[ 'Minimum Light Intensity' ].layout.display = None
            Widgets[ 'Maximum Light Intensity' ].layout.display = None
            Widgets[ 'Number of Light Intensities' ].layout.display = None
            Widgets[ 'Light Intensity Sample Type'].layout.display = None
    
                
        if GUI_Mode == 'Varied Light Source':
            
            Widgets[ 'Light Intensity' ].layout.display = None
    
        if GUI_Mode == 'Typical Indoor Day':
                
            Widgets[ 'Indoor Day Start Time (Hour)' ].layout.display = None
            Widgets[ 'Indoor Day End Time (Hour)' ].layout.display = None
            Widgets[ 'Indoor Day Peak Intensity' ].layout.display = None
            
        if GUI_Mode == 'Many-System Comparison':
            Widgets[ "Junction Type" ].value = 'Single-Junction' # Don't allow MJ in this case
            Widgets[ "Junction Type" ].layout.display = 'none'
            Widgets[ 'Non-Radiative Loss Mode' ].layout.display = 'none'
            Widgets[ "Material Model" ].layout.display = 'none'
            Widgets[ "Simulated EQE Model" ].layout.display = 'none'
            Widgets[ 'E_g (eV)' ].layout.display = 'none'
            Widgets[ 'Above Gap EQE' ].layout.display = 'none'
            Widgets[ 'Below Gap EQE' ].layout.display = 'none'
            Widgets[ 'Urbach Energy (eV)' ].layout.display = 'none'
            Widgets[ 'Excitonic Disorder (eV)' ].layout.display = 'none'
            Widgets[ "Experimental EQE Model" ].layout.display = 'none'
    
        if GUI_Mode == "Multi-Junction - One Fixed Cell":
            
            Widgets[ "Junction Type" ].value = 'Multi-Junction' # Ensure MJ in this case
            Widgets[ "Junction Type" ].layout.display = 'none'
            Widgets[ 'Fixed Multi-Junction Cell' ].layout.display = None
            
            if Widgets[ 'Fixed Multi-Junction Cell' ].value == 'Top':
                
                # Need to hide the appropriate widgets to allow bottom cell simulation
                
                if Widgets[ "Bottom Cell Material Model" ].value == 'Simulated':
                    
                    Gap_Inputs_Label.layout.display = None
                    Widgets[ 'Bottom Cell E_g (eV)' ].layout.display = 'none'     
                    Widgets[ "Bottom Cell Experimental EQE Model" ].layout.display = 'none'
                    
                    Widgets[ 'Minimum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Maximum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Number of E_gs' ].layout.display = None
                    
                    Widgets[ "Bottom Cell Simulated EQE Model" ].layout.display = None
                    Widgets[ 'Bottom Cell Above Gap EQE' ].layout.display = None
                    
                    if Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Step Function':
                        
                        Widgets[ 'Bottom Cell Below Gap EQE' ].layout.display = None
    
                    elif Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Urbach Tail':
                        
                        Widgets[ 'Bottom Cell Urbach Energy (eV)' ].layout.display = None
                        
                    elif Widgets[ "Bottom Cell Simulated EQE Model" ].value == 'Disordered Organic PV':
                        
                        Widgets[ 'Bottom Cell Excitonic Disorder (eV)' ].layout.display = None
                            
                else:
                    
                    Widgets[ "Bottom Cell Experimental EQE Model" ].layout.display = 'none'  # Use all EQEs
                
            if Widgets[ 'Fixed Multi-Junction Cell' ].value == 'Bottom':
                
                # Need to hide the appropriate widgets to allow top cell simulation
                
                if Widgets[ "Material Model" ].value == 'Simulated':
                    
                    Gap_Inputs_Label.layout.display = None
                    Widgets[ 'E_g (eV)' ].layout.display = 'none'     
                    Widgets[ "Experimental EQE Model" ].layout.display = 'none'
                    
                    Widgets[ 'Minimum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Maximum E_g (eV)' ].layout.display = None        
                    Widgets[ 'Number of E_gs' ].layout.display = None
                    
                    Widgets[ "Simulated EQE Model" ].layout.display = None
                    Widgets[ 'Above Gap EQE' ].layout.display = None
                    
                    if Widgets[ "Simulated EQE Model" ].value == 'Step Function':
                        
                        Widgets[ 'Below Gap EQE' ].layout.display = None
    
                    elif Widgets[ "Simulated EQE Model" ].value == 'Urbach Tail':
                        
                        Widgets[ 'Urbach Energy (eV)' ].layout.display = None
                        
                    elif Widgets[ "Simulated EQE Model" ].value == 'Disordered Organic PV':
                        
                        Widgets[ 'Excitonic Disorder (eV)' ].layout.display = None
                            
                else:
                    
                    Widgets[ "Experimental EQE Model" ].layout.display = 'none' 
                    
          #          pass # Widgets[ "Bottom Cell Experimental EQE Model" ].layout.display = None # Use all EQEs
                
                
    for Key in list( Widgets.keys() ):
        
        Widgets[ Key ].observe( GUI_Updater, names = 'value' )
    
    # NEED FUNCTION FOR CHANGING MULTI-TO-SINGLE JUNCTION
    # NEED FUNCTION FOR CHANGING EQE SIMULAITON MODE
    
    #==========================================================================
    # Compile GUI
    #==========================================================================
    
    Gap = Label( '' )
    Top_Cell_Label = Label( '𝘛𝘰𝘱 𝘊𝘦𝘭𝘭 𝘐𝘯𝘱𝘶𝘵𝘴:' )
    Bottom_Cell_Label = Label( '𝘉𝘰𝘵𝘵𝘰𝘮 𝘊𝘦𝘭𝘭 𝘐𝘯𝘱𝘶𝘵𝘴:' )
    Gap_Inputs_Label = Label( '𝘎𝘢𝘱 𝘐𝘯𝘱𝘶𝘵𝘴:' )
    Top_Cell_Label.layout.display = 'none'
    Bottom_Cell_Label.layout.display = 'none'
    Gap_Inputs_Label.layout.display = 'none'
    
    Layer_Inputs = []
    
    for j in range( Widgets[ "Number of Layers in Device" ].max ):
        
        Layer_Inputs.append( Widgets[ "Layer Type " + str( j + 1 ) ] )
        Layer_Inputs.append( Widgets[ "Layer Thickness " + str( j + 1 ) + " (nm)" ] )
    
    Inputs_Box = VBox( [ 
        
        Label( value = '𝗖𝗼𝗻𝘁𝗿𝗼𝗹 𝗣𝗮𝗿𝗮𝗺𝗲𝘁𝗲𝗿𝘀' ),
        
        Widgets[ 'Simulation Mode' ],      
        Gap,         
        Widgets[ 'Optically Modelled' ],
        Widgets[ 'Junction Type' ],
        Widgets[ 'Fixed Multi-Junction Cell' ],     
        Widgets[ 'Multi-Junction Terminal Type' ],          
        Widgets[ 'Radiative Coupling Top to Bottom' ],
        Widgets[ 'Radiative Coupling Bottom to Top' ],
        Widgets[ 'Module Simulation' ],     
        Widgets[ 'Position Spacing (nm)' ],
        Widgets[ 'Generate Plots' ],
        Gap,          
        
        #----------------------------------------------------------------------
        Label( value = '𝗠𝗮𝘁𝗲𝗿𝗶𝗮𝗹 𝗣𝗮𝗿𝗮𝗺𝗲𝘁𝗲𝗿𝘀' ),
        #----------------------------------------------------------------------
        Widgets[ "Surrounding Medium" ],
        Widgets[ "Incident Angle (deg)" ],
        Widgets[ "Number of Layers in Device" ],
        Widgets[ "Incoherent Substrate" ],
        Widgets[ "Substrate Index" ],
        Widgets[ "Infinite Final Layer" ],
        Widgets[ "Active Layer Index" ],
        Widgets[ "Bottom Active Layer Index" ],
        Widgets[ "Internal Quantum Efficiency" ],
        Widgets[ "Bottom Cell Internal Quantum Efficiency" ]]
        
        + Layer_Inputs
    
        + [ Top_Cell_Label,     
        Widgets[ 'Material Model' ],
        Widgets[ 'Simulated EQE Model' ],
        Widgets[ 'Experimental EQE Model' ],
        Widgets[ 'E_g (eV)' ],
        Widgets[ 'Above Gap EQE' ],
        Widgets[ 'Below Gap EQE' ],
        Widgets[ 'Urbach Energy (eV)' ],
        Widgets[ 'Excitonic Disorder (eV)' ],
        Widgets[ 'Non-Radiative Loss Mode' ],
        Bottom_Cell_Label, 
        Widgets[ 'Bottom Cell Material Model' ],
        Widgets[ 'Bottom Cell Simulated EQE Model' ],
        Widgets[ 'Bottom Cell E_g (eV)' ],
        Widgets[ 'Bottom Cell Above Gap EQE' ],
        Widgets[ 'Bottom Cell Below Gap EQE' ],
        Widgets[ 'Bottom Cell Urbach Energy (eV)' ],
        Widgets[ 'Bottom Cell Excitonic Disorder (eV)' ],    
        Widgets[ "Bottom Cell Experimental EQE Model" ],
        Widgets[ 'Bottom Cell Non-Radiative Loss Mode' ],
        Gap_Inputs_Label,
        Widgets[ 'Minimum E_g (eV)' ],
        Widgets[ 'Maximum E_g (eV)' ],
        Widgets[ 'Number of E_gs' ],
        Gap,
        
        #----------------------------------------------------------------------
        Label( value = '𝗗𝗲𝘃𝗶𝗰𝗲 𝗣𝗮𝗿𝗮𝗺𝗲𝘁𝗲𝗿𝘀' ),
        #----------------------------------------------------------------------
        
        Widgets[ 'Number of Cells' ],
        Widgets[ 'Module Width (cm)' ],
        Widgets[ 'Module Length (cm)' ],
        Widgets[ 'Area (cm2)' ],
        Widgets[ 'Interconnect Width (cm)' ],
        Widgets[ 'R_series (Ohms)' ],
        Widgets[ 'R_shunt (Ohms cm2)' ],
        Widgets[ 'Bottom Cell R_series (Ohms)' ],
        Widgets[ 'Bottom Cell R_shunt (Ohms cm2)' ],
        Widgets[ "Include Point Defects" ],
        Widgets[ "Simulate Number of Defects" ],
        Widgets[ "Number of Defects" ],
        Widgets[ 'Point Defect Resistance (Ohms cm2)' ],
        Widgets[ 'Point Defect Proclivity f_x (1/cm)' ], 
        Widgets[ 'Point Defect Proclivity f_y (1/cm)' ],
        Widgets[ 'A_onset (cm)' ],
    
        Gap,
        #----------------------------------------------------------------------
        Label( value = '𝗟𝗶𝗴𝗵𝘁 𝗣𝗮𝗿𝗮𝗺𝗲𝘁𝗲𝗿𝘀' ),
        #----------------------------------------------------------------------
    
        Widgets[ 'Light Source' ],
        Widgets[ 'Light Intensity' ],
        Widgets[ 'Minimum Light Intensity' ],
        Widgets[ 'Maximum Light Intensity' ],
        Widgets[ 'Number of Light Intensities' ],
        Widgets[ 'Light Intensity Sample Type'],
        Widgets[ 'Indoor Day Start Time (Hour)' ],
        Widgets[ 'Indoor Day End Time (Hour)' ],
        Widgets[ 'Indoor Day Peak Intensity' ],
        Widgets[ 'Intensity Unit' ],
        Widgets[ 'Minimum Energy (eV)' ],
        Widgets[ 'Maximum Energy (eV)' ],
        Widgets[ 'Number of Energies' ],
        Gap,
        
        #----------------------------------------------------------------------
        Label( value = '𝗔𝗺𝗯𝗶𝗲𝗻𝘁 𝗣𝗮𝗿𝗮𝗺𝗲𝘁𝗲𝗿𝘀' ),
        #----------------------------------------------------------------------
        Widgets[ 'Temperature' ],
        Widgets[ 'Temperature Unit' ],
        
        
    ] )
    
    Intialise_Button = Button()
    
    Intialise_Button.on_click( GUI_Updater )
    
    Intialise_Button.click()
    
    #-------------------------------------------------------------------------
    # Create Output Box
    #-------------------------------------------------------------------------
    
    Defaults_Button = Button( description = 'Reset Default Values' )
    Run_Simulations_Button = Button( description = 'Run Simulations' )
    Copy_Output = Button( description = 'Copy Output' )
    
    def Default_Value_Loader( Button ):
        
        Widgets[ "Light Source" ].value = 'AM1.5G'
        
        for j in range( len( Dropdown_Widget_DF ) ):
            Row = Dropdown_Widget_DF.iloc[ j ]
            Widgets[ Row[ "Description" ] ].value = str( Row[ "Default Value" ] )
            
        for j in range( len( Float_Widget_DF ) ):
            Row = Float_Widget_DF.iloc[ j ]
            Widgets[ Row[ "Parameter" ] ].value = Row[ "Default Value" ]
            
        for j in range( len( Int_Widget_DF ) ):
            Row = Int_Widget_DF.iloc[ j ]
            Widgets[ Row[ "Parameter" ] ].value = Row[ "Default Value" ]
        
        Value_Storer( '' )
    
    def Value_Storer( Button ):
        
    #    Widgets[ "Light Source" ].value = 'AM1.5G'
        
        for j in range( len( Dropdown_Widget_DF ) ):
            Row = Dropdown_Widget_DF.iloc[ j ]
            Row[ "Saved Value" ] = Widgets[ Row[ "Description" ] ].value
            Dropdown_Widget_DF.iloc[ j ] = Row
            
        for j in range( len( Float_Widget_DF ) ):
            Row = Float_Widget_DF.iloc[ j ]
            Row[ "Saved Value" ] = Widgets[ Row[ "Parameter" ] ].value
            Float_Widget_DF.iloc[ j ] = Row       
            
        for j in range( len( Int_Widget_DF ) ):
            Row = Int_Widget_DF.iloc[ j ]
            Row[ "Saved Value" ] = Widgets[ Row[ "Parameter" ] ].value
            Int_Widget_DF.iloc[ j ] = Row  
            
        Excel_Writer = ExcelWriter( Parameter_Filepath, engine = 'xlsxwriter' )   
        Excel_Workbook = Excel_Writer.book
        Dropdown_Sheet = Excel_Workbook.add_worksheet( 'Dropdowns' ) 
        Float_Sheet = Excel_Workbook.add_worksheet( 'Floats' ) 
        Int_Sheet = Excel_Workbook.add_worksheet( 'Integers' ) 
    
        Dropdown_Widget_DF.to_excel( Excel_Writer, sheet_name = 'Dropdowns', index = False )   
        Float_Widget_DF.to_excel( Excel_Writer, sheet_name = 'Floats', index = False )   
        Int_Widget_DF.to_excel( Excel_Writer, sheet_name = 'Integers', index = False )   
        
        Dropdown_Sheet.autofit()
        Float_Sheet.autofit()
        Int_Sheet.autofit()
        
        Excel_Workbook.close()
        
    def GUI_Graph_Generator( Change ):
        
        #--------------------------------------------------------------------------
        # Determine Available Graphs
        #--------------------------------------------------------------------------
        
        Output_Path = path.join( Current_Working_Directory, "Output" )
    
        Available_Outputs = listdir( Output_Path )
        
        Recent_Output_Path = path.join( Output_Path, Available_Outputs[ -1 ] )
        
        global Available_Graphs 
        
        Available_Graphs = listdir( Recent_Output_Path )
        
        Available_Graphs = [ Graph.replace( '.jpg', '' ) for Graph in Available_Graphs if '.jpg' in Graph ]
    
        #--------------------------------------------------------------------------
        # Generate and Compile Widgets
        #--------------------------------------------------------------------------
        
        global Graph_Selector 
        
        Graph_Selector = RadioButtons( options = Available_Graphs,
                                      
                                       value = Available_Graphs[ 0 ] )
    
        global Images
        
        Images = { Graph : Output_Widget() for Graph in Available_Graphs }
        
        for Graph in Available_Graphs:
            
            with Images[ Graph ]:
                
                display( Image( path.join( Recent_Output_Path, Graph + '.jpg' ) ,
                        
                        width = 600 ) )
                        
        Graph_Selector.observe( GUI_Graph_Updater, names = 'value' )
        
        if 'PCE' in Available_Graphs:

            Graph_Selector.value = 'PCE' 
        
        else:
            
            Graph_Selector.value = 'PhotonFluxes' 
            
        Graph_Box = VBox( [ Graph_Selector ] + [ Images[ Graph ]
                                                
                                                for Graph in Available_Graphs ] )
        
        #--------------------------------------------------------------------------
        # Update GUI
        #--------------------------------------------------------------------------
            
        Output_Box.children = Output_Box.children[ :-1 ] + ( Graph_Box, )
        
    def GUI_Graph_Updater( Change ):
        
        New_Graph = Change[ 'new' ]
        
        for Graph in Available_Graphs:
        
            if Graph == New_Graph:
                
                Images[ Graph ].layout.display = None
            
            else:
                
                Images[ Graph ].layout.display = "none"
                
    def Output_Copier( Button ):
    
        DataFrame( Output ).to_clipboard( index = False )
    
    Copy_Output.on_click( Output_Copier )
    
    Copy_Output.disabled = True
    
    #--------------------------------------------------------------------------                
    # Compile GUI
    #--------------------------------------------------------------------------
    
    Defaults_Button.on_click( Default_Value_Loader )
    
    Buttons_Box = HBox( [ Defaults_Button, Run_Simulations_Button, Copy_Output ] )
    
    Output_Box = VBox( [ Buttons_Box, Gap ])
    
    User_Interface = HBox( [ Inputs_Box, Gap, Output_Box ] ) 
    
#%%==========================================================================10
# 10. Run the Simulations
#==============================================================================

"""

Note, for multi-junction simulations with a varied top/bottom cell, the inputs 

are assumed to be the inputs of the fixed cell. Whereas in the case of 'Single

Simulation', 'Varied Intensity', 'Varied Light Source', or 'Typical Indoor Day', 

they are assumed to be the inputs of the top cell.

"""
    
def Simulation_Controller():
        
    #==============================================================================
    # Part 0 - Create Output File Using Data and Time
    #==============================================================================
    
    global Output 
    
    Output_Path = path.join( Current_Working_Directory, "Output" )
    
    if 'Output' not in listdir( Current_Working_Directory ):
        
        mkdir( Output_Path )
    
    Date, Time = Date_Time_Collector() 

    Output_Path = path.join( Output_Path, "Output_" + Date + '_' + Time  )

    mkdir( Output_Path )

    #==============================================================================
    # Part i - Load Tempearture and Photon Fluxes
    #==============================================================================
    
    Parameters_Path = path.join( Current_Working_Directory, "Parameter_Values.xlsx" )
    
    Dropdown_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Dropdowns', engine='openpyxl' )
    Float_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Floats', engine='openpyxl')
    Int_Widget_DF = read_excel( Parameter_Filepath, sheet_name = 'Integers', engine='openpyxl' )
    
    global Parameters
    Parameters = {}
        
    for j in range( len( Dropdown_Widget_DF ) ):
        Row = Dropdown_Widget_DF.iloc[ j ]
        Parameters[ Row[ "Description" ] ] = Row[ "Saved Value" ] 
        
    for j in range( len( Float_Widget_DF ) ):
        Row = Float_Widget_DF.iloc[ j ]
        Parameters[ Row[ "Parameter" ] ] = Row[ "Saved Value" ] 
        
    for j in range( len( Int_Widget_DF ) ):
        Row = Int_Widget_DF.iloc[ j ]
        Parameters[ Row[ "Parameter" ] ] = Row[ "Saved Value" ] 
    
    Parameter_Storing_Path = path.join( Output_Path, "Parameters.json" )
    
    with open( Parameter_Storing_Path, 'w' ) as File:
        json.dump( Parameters, File, indent = 2, cls = NpEncoder )
    File.close()
                                                                      
    #------------------------------------------------------------------------------
    # Energies and Wavelengths:
    #------------------------------------------------------------------------------
    
    min_Energy = Parameters[ "Minimum Energy (eV)" ]
    max_Energy = Parameters[ "Maximum Energy (eV)" ]
    N_Energy = Parameters[ "Number of Energies" ]
    
    Energies = linspace( min_Energy, max_Energy, N_Energy )
    Energies = Energies[::-1]  # invert order to get increasing wavelengths
    Wavelengths = Energy_Wavelength_Converter( Energies )
    
    #------------------------------------------------------------------------------
    # Temperature:
    #------------------------------------------------------------------------------
    
    T = Parameters[ "Temperature" ]
    Temperature_Unit = Parameters[ "Temperature Unit" ] 
    
    if Temperature_Unit == 'C':
        T += 273.15
    if Temperature_Unit == 'F':
        T = ( T - 32 ) * 5 / 9 + 273.15
        
    kT = k * T / e
    
    Black_Body_Fluxes = Planck_Photon_Flux_Wavelength( Wavelengths, kT )  # Photons/m2/s/nm
    
    One_Sun_Photon_Fluxes = Photon_Flux_Interpolator_and_Scaler( "AM1.5G",
                                                                 Wavelengths,
                                                                 100.3 ) # Photons /m2/s/nm
    
    #==============================================================================
    # Load Photon Fluxes (Single Spectrum, Single Intensity)
    #==============================================================================
    
    Mode = Parameters[ "Simulation Mode" ]
    
    if Mode not in Simulation_Modes:
        
        print( "Simulation mode not valid, please update file." )
        
        exit()
        
    Generate_Plots = Parameters[ "Generate Plots" ]
    
    if Generate_Plots == 'True':
        Generate_Plots = True
    if Generate_Plots == 'False':
        Generate_Plots = False
        
    if Mode in [ "Single Simulation", "Many-System Comparison", "Varied Gap", 
                
                "Varied Lower Limit", "Multi-Junction - One Fixed Cell" ]:
        
        Light_Source = Parameters[ "Light Source" ] 
        
        Light_Power = Parameters[ "Light Intensity" ]
        
        Intensity_Unit = Parameters[ "Intensity Unit" ]
    
        if Intensity_Unit == 'W/m2':    
            Light_Power /= 10 
        if Intensity_Unit == 'lux':
            Lux_Constant = Lux_Constant_Calculator( Light_Source, Wavelengths, Energies )
            Light_Power = Light_Power / Lux_Constant / 10 # mW/cm2
    
        Incident_Photon_Fluxes = Photon_Flux_Interpolator_and_Scaler( Light_Source,
                                                                     Wavelengths,
                                                                     Light_Power ) # Photons /m2/s/nm
    
        if Generate_Plots:
            
            plt.figure( dpi = 200 )
            plt.plot( Wavelengths, Incident_Photon_Fluxes )
            plt.xlabel( 'Wavelength, (nm)' )
            plt.ylabel( 'Incident Photon Flux (Photons m$^{-2}$s$^{-1}$nm$^{-1}$)' )
            plt.xlim( [ 300, 4000 ] );
            plt.savefig( path.join( Output_Path, 'PhotonFluxes.jpg' ) )
            
    #------------------------------------------------------------------------------
    # Load Photon Fluxes (Single Spectrum, Varied Intensity)
    #------------------------------------------------------------------------------
    
    elif Mode == 'Varied Intensity':
        
        Light_Source = Parameters[ "Light Source" ]     
        Intensity_Unit = Parameters[ "Intensity Unit" ]
        
        Min_Light_Power = Parameters[ "Minimum Light Intensity" ]
        Max_Light_Power = Parameters[ "Maximum Light Intensity" ]
        N_Light_Power = Parameters[ "Number of Light Intensities" ]
        Light_Power_Sample_Type = Parameters[ "Light Intensity Sample Type" ]
            
        if Light_Power_Sample_Type == 'Linear':
            
            Light_Powers_initial = linspace( Min_Light_Power, Max_Light_Power, N_Light_Power )
    
        if Light_Power_Sample_Type == 'Logarithmic':
                
            Light_Powers_initial = logspace( log10( Min_Light_Power ), log10( Max_Light_Power ), N_Light_Power )
    
        if Intensity_Unit == 'mW/cm2':
            Light_Powers = Light_Powers_initial.copy()
            
        if Intensity_Unit == 'W/m2':    
            Light_Powers = Light_Powers_initial.copy() / 10 
            
        if Intensity_Unit == 'lux':
            Lux_Constant = Lux_Constant_Calculator( Light_Source, Wavelengths, Energies )
            Light_Powers = Light_Powers_initial.copy() / Lux_Constant / 10 # mW/cm2
    
        Incident_Photon_Fluxes = { Light_Power : Photon_Flux_Interpolator_and_Scaler( Light_Source,
                                                                     Wavelengths,
                                                                     Light_Power ) # Photons /m2/s/nm
                                  for Light_Power in Light_Powers }
        
    #------------------------------------------------------------------------------
    # Load Photon Fluxes (Varied Spectrum, Single Intensity)
    #------------------------------------------------------------------------------
    
    elif Mode == 'Varied Light Source':
                    
        Light_Power = Parameters[ "Light Intensity" ]
        
        Intensity_Unit = Parameters[ "Intensity Unit" ]
    
        Select_Light_Sources = [ Light_Source for Light_Source in Light_Sources if 'Time' not in Light_Source ]

        Light_Powers_initial = { Light_Source : Light_Power for Light_Source in Select_Light_Sources }

        if Intensity_Unit == 'mW/cm2':
            Light_Powers = { Light_Source : Light_Power for Light_Source in Select_Light_Sources }
        if Intensity_Unit == 'W/m2':    
            Light_Powers = { Light_Source : Light_Power / 10 for Light_Source in Select_Light_Sources }
        if Intensity_Unit == 'lux':
            # Fixed Illuminance
            Light_Powers = { Light_Source : Light_Power / Lux_Constant_Calculator( Light_Source, Wavelengths, Energies ) / 10                        
                            for Light_Source in Select_Light_Sources }
    
        Incident_Photon_Fluxes = { Light_Source : Photon_Flux_Interpolator_and_Scaler( Light_Source,
                                                                     Wavelengths,
                                                                     Light_Powers[ Light_Source ] )
                                                                     
                                                                     for Light_Source in Select_Light_Sources } # Photons /m2/s/nm
        
        if Generate_Plots:
            
            plt.figure( dpi = 200 )
            for Light_Source in Select_Light_Sources:    
                plt.plot( Wavelengths, Incident_Photon_Fluxes[ Light_Source ], label = Light_Source )
            plt.xlabel( 'Wavelength, $\\lambda$ (nm)' )
            plt.ylabel( 'Incident Photon Flux (Photons m$^{-2}$s$^{-1}$nm$^{-1}$)' )
            plt.xlim( [ 300, 4000 ] );
            plt.savefig( path.join( Output_Path, 'PhotonFluxes.jpg' ) )

    elif Mode == 'Typical Indoor Day':
        
        # Specify Time Range
        
        Start_Time = Parameters[ "Indoor Day Start Time (Hour)" ]
        End_Time = Parameters[ "Indoor Day End Time (Hour)" ]
    
        Day_Light_Sources = [ 'Time_' + str( j ) + '00' if len( str( j ) ) == 2 
                             
                             else 'Time_0' + str( j ) + '00' 
                             
                             for j in range( int( Start_Time ), int( End_Time + 1 ) ) ]
             
        #==========================================================================
        # Determine UNSCALED Intensity for all spectra
        #==========================================================================
    
        Unscaled_Intensities = {}   
        
        for Light_Source in Day_Light_Sources: 
    
            Wavelengths_j, Fluxes_j = Photon_Flux_Spectrum_Loader( Light_Source )
            
            Energies_j = Energy_Wavelength_Converter( Wavelengths_j )  # in eV
            
            P_light = Light_Power_Calculator(Wavelengths_j,
            
                                             Energies_j,
            
                                             Fluxes_j)
    
            print( Light_Source, ': Unscaled P_light = ', P_light, 'mW/cm2' )
    
            Unscaled_Intensities[ Light_Source ] = P_light
        
        #==========================================================================
        # Determine Scale Factor
        #==========================================================================
        
        Peak_Intensity_Index = where( list( Unscaled_Intensities.values() ) == max( list( Unscaled_Intensities.values() ) ) )[ 0 ][ 0 ]
        Peak_Intensity = list( Unscaled_Intensities.values() )[ Peak_Intensity_Index ]
    
        Peak_Time_Light_Source = list( Unscaled_Intensities.keys() )[ Peak_Intensity_Index ]
    
        Desired_Peak_Intensity = Parameters[ "Indoor Day Peak Intensity" ]
        Intensity_Unit = Parameters[ "Intensity Unit" ]
        
        if Intensity_Unit == 'mW/cm2':
            Desired_Peak_Intensity = Desired_Peak_Intensity
            
        if Intensity_Unit == 'W/m2':    
            Desired_Peak_Intensity = Desired_Peak_Intensity / 10 
            
        if Intensity_Unit == 'lux':
            Constant_of_Prop = Lux_Constant_Calculator( Peak_Time_Light_Source, Wavelengths, Energies )
            Desired_Peak_Intensity = Desired_Peak_Intensity/ Constant_of_Prop / 10 # mW/cm2
    
        Scale_Factor = Peak_Intensity / Desired_Peak_Intensity # Unitless, need to divide Light_Powers by this
    
        #==========================================================================
        # Determine Photon Fluxes
        #==========================================================================
        
        Light_Powers = { Light_Source : Unscaled_Intensities[ Light_Source ] / Scale_Factor 
                        
                        for Light_Source in Day_Light_Sources } # mW/cm2
        
        Incident_Photon_Fluxes = { Light_Source : Photon_Flux_Interpolator_and_Scaler( Light_Source,
                                                                     Wavelengths,
                                                                     Light_Powers[ Light_Source ] ) # Photons /m2/s/nm
                                  
                                  for Light_Source in Day_Light_Sources }    
    
        if Generate_Plots: 
            
            plt.figure( dpi = 200 )
            
            for Light_Source in Day_Light_Sources:
                
                plt.plot( Wavelengths, Incident_Photon_Fluxes[ Light_Source ], label = Light_Source )
            plt.xlabel( 'Wavelength (nm)' )
            plt.ylabel( 'Photons/s/m2/nm' )
            plt.legend()
            plt.xlim( [ 200, 2000 ] )
            plt.savefig( path.join( Output_Path, 'PhotonFluxes.jpg' ) )                
            
    #==============================================================================
    # Part ii - Load Device and Material Parameters
    #==============================================================================
    
    #------------------------------------------------------------------------------
    # Device Parameters:
    #------------------------------------------------------------------------------
    
    Junction_Type = Parameters[ "Junction Type" ]
    
    if Junction_Type == 'Single-Junction':
        
        Multijunction_Simulation = False
        
    else:
        Multijunction_Simulation = True 
    
    Module_Simulation = Parameters[ "Module Simulation" ]
    
    #--------------------------------------------
    # Model Areas:
    #--------------------------------------------
    
    if Module_Simulation == 'False': 
        
        Area_cell = Parameters[ "Area (cm2)" ]
    
    else:
            
        L_x = Parameters[ "Module Width (cm)" ] # Interconnects are assumed to run along the y-plane, with delta along x 
        L_y = Parameters[ "Module Length (cm)" ] 
        Area_Module = L_x * L_y # cm2
        delta = Parameters[ "Interconnect Width (cm)" ] # cm
        N_Cells = Parameters[ "Number of Cells" ]
        
        if N_Cells  < 1 + L_x / delta:
                
            Area_cell = ( Area_Module - ( N_Cells - 1  ) * delta * L_y ) / N_Cells 
            
        else:
            
            Area_cell = 0 
        
    #--------------------------------------------    
    # Model Parasitic Resistances:
    #--------------------------------------------
        
    R_s = Parameters[ "R_series (Ohms)" ]
    R_sh = Parameters[ "R_shunt (Ohms cm2)" ]
    
    if R_sh == 'inf':
        R_sh = inf
    
    else:
        R_sh = float( R_sh )

    #--------------------------------------------    
    # MULTI-JUNCTION: Model Parasitic Resistances:
    #--------------------------------------------
        
    if Multijunction_Simulation:
        
        BC_R_s = Parameters[ "Bottom Cell R_series (Ohms)" ]
        BC_R_sh = Parameters[ "Bottom Cell R_shunt (Ohms cm2)" ]
    
        if BC_R_sh == 'inf':
            BC_R_sh = inf
        
        else:
            BC_R_sh = float( BC_R_sh )    
        
    #---------------------------------------------
    # Account for Point Defects
    #--------------------------------------------
        
    if Parameters[ "Include Point Defects" ] == 'True':
        
        r_Defect = Parameters[ "Point Defect Resistance (Ohms cm2)" ]
        
        if Parameters[ "Simulate Number of Defects" ] == 'False':
            
            N_Defects = Parameters[ "Number of Defects" ]
            
        else:
            
            f_x, f_y = Parameters[ 'Point Defect Proclivity f_x (1/cm)' ], Parameters[ 'Point Defect Proclivity f_y (1/cm)' ]
    
            A_onset = Parameters[ 'A_onset (cm)' ]
            
            N_Defects = int( floor( L_x * L_y / A_onset * exp( f_x * L_x + f_x * L_y ) ) )
    
        ### Distribute defects if using a module:
            
        if Module_Simulation == 'True':
            
            N_Defects = int( N_Defects / N_Cells ) # per cell
    
        R_sh = 1 / ( 1 / R_sh + N_Defects / r_Defect )  # effective shunt/parallel resistance
        
    #------------------------------------------------------------------------------
    # Material Parameters (Single EQE)
    #------------------------------------------------------------------------------
    
    if Mode in [ "Single Simulation", 'Varied Intensity', 'Varied Light Source', 
                
                'Typical Indoor Day' ]:
    
        if Parameters[ "Optically Modelled" ] != "True":
            
            EQE_Mode = Parameters[ "Material Model" ]
                    
            #--------------------------------------------------------------------------
            if EQE_Mode == "Experimental (Loaded)":
            #--------------------------------------------------------------------------
                
                EQE_Filename = Parameters[ "Experimental EQE Model" ]
                
                if '.xlsx' in EQE_Filename:
                    
                    EQE_Filename.replace( '.xlsx', '')
                    
                E_g = Systems_Database[ "E_g [eV]" ][ EQE_Filename ]
                E_lower = Systems_Database[ "E_lower [eV]" ][ EQE_Filename ]
                V_oc_exp = Systems_Database[ 'V_oc [V]' ][ EQE_Filename ]
        
                EQEs = EQE_Loader_and_Interpolator( EQE_Filename + '.xlsx', Energies, E_lower )
            
                # print( EQE_Filename, ': E_g = ', E_g, 'eV; V_oc = ', V_oc_exp, 'V' )            
            
            #--------------------------------------------------------------------------
            if EQE_Mode == "Simulated":
            #--------------------------------------------------------------------------
                    
                Simulated_EQE_Type = Parameters[ "Simulated EQE Model" ]
            
                E_g = Parameters[ "E_g (eV)" ]
            
                EQE_max = Parameters[ "Above Gap EQE" ]
                
                if Simulated_EQE_Type == "Step Function":
                
                    EQE_min = Parameters[ "Below Gap EQE" ]
            
                    EQEs = Step_Function_EQE_Simulator( Energies, E_g, EQE_min, EQE_max )
            
                if Simulated_EQE_Type == "Urbach Tail":
                    
                    E_U = Parameters[ "Urbach Energy (eV)" ]
        
                    EQEs = Urbach_Tail_EQE_Simulator( Energies, E_g, E_U, EQE_max )
            
                if Simulated_EQE_Type == "Disordered Organic PV":
                    
                    sigma_s = Parameters[ "Excitonic Disorder (eV)" ]
        
                    EQEs = Exction_Absorption_EQE_Simulator( Energies, E_g, sigma_s, EQE_max, kT )
                 
            if Multijunction_Simulation:
            
                #------------------------------------------------------------------
                # Assume Above EQEs are for Top Cell; Load EQEs for Bottom Cell
                #------------------------------------------------------------------            
                
                Bottom_Cell_EQE_Mode = Parameters[ "Bottom Cell Material Model" ]
            
                #------------------------------------------------------------------
                if Bottom_Cell_EQE_Mode == "Experimental (Loaded)":
                #------------------------------------------------------------------
                    
                    BC_EQE_Filename = Parameters[ "Bottom Cell Experimental EQE Model" ]
                    
                    if '.xlsx' in BC_EQE_Filename:
                        
                        BC_EQE_Filename.replace( '.xlsx', '')
                        
                    BC_E_g = Systems_Database[ "E_g [eV]" ][ BC_EQE_Filename ]
                    BC_E_lower = Systems_Database[ "E_lower [eV]" ][ BC_EQE_Filename ]
                    BC_V_oc_exp = Systems_Database[ 'V_oc [V]' ][ BC_EQE_Filename ]
            
                    BC_EQEs = EQE_Loader_and_Interpolator( BC_EQE_Filename + '.xlsx', Energies, BC_E_lower )
                
                    print( 'Bottom Cell; ', BC_EQE_Filename, ': E_g = ', BC_E_g, 'eV; V_oc = ', BC_V_oc_exp, 'V' )            
                
                #------------------------------------------------------------------
                if Bottom_Cell_EQE_Mode == "Simulated":
                #------------------------------------------------------------------
                        
                    BC_Simulated_EQE_Type = Parameters[ "Bottom Cell Simulated EQE Model" ]
                
                    BC_E_g = Parameters[ "Bottom Cell E_g (eV)" ]
                
                    BC_EQE_max = Parameters[ "Bottom Cell Above Gap EQE" ]
                    
                    if BC_Simulated_EQE_Type == "Step Function":
                    
                        BC_EQE_min = Parameters[ "Bottom Cell Below Gap EQE" ]
                
                        BC_EQEs = Step_Function_EQE_Simulator( Energies, BC_E_g, BC_EQE_min, BC_EQE_max )
                
                    if BC_Simulated_EQE_Type == "Urbach Tail":
                        
                        BC_E_U = Parameters[ "Bottom Cell Urbach Energy (eV)" ]
            
                        BC_EQEs = Urbach_Tail_EQE_Simulator( Energies, BC_E_g, BC_E_U, BC_EQE_max )
                
                    if BC_Simulated_EQE_Type == "Disordered Organic PV":
                        
                        BC_sigma_s = Parameters[ "Bottom Cell Excitonic Disorder (eV)" ]
            
                        BC_EQEs = Exction_Absorption_EQE_Simulator( Energies, BC_E_g, BC_sigma_s, BC_EQE_max, kT )
                    
                if Generate_Plots:
                    plt.figure( dpi = 200 )
                    plt.plot( Energies, EQEs, label = 'Top Cell' )
                    if Multijunction_Simulation:
                        plt.plot( Energies, BC_EQEs, label = 'Bottom Cell' )
                        plt.legend()
                        
                    plt.xlabel( '$E$ (eV)' )
                    plt.ylabel( '$\\mathrm{EQE}_{\\mathrm{PV}}$' )
                    plt.savefig( path.join( Output_Path, 'EQE.jpg' ) )
                
        else:
                                
            Surrounding_Medium = Parameters[ "Surrounding Medium" ]
            Incident_Angle = Parameters[ "Incident Angle (deg)" ]

            Number_of_Layers_in_Device = Parameters[ "Number of Layers in Device" ]
            
            Layer_Types = [ Parameters[ "Layer Type " + str( j + 1 ) ] for j in range( Number_of_Layers_in_Device ) ]              
            Layer_Thicknesses = [ Parameters[ "Layer Thickness " + str( j + 1 ) + " (nm)" ] for j in range( Number_of_Layers_in_Device ) ]              
            Active_Layer_Index = Parameters[ "Active Layer Index" ]
            IQEs = array( [ Parameters[ "Internal Quantum Efficiency" ] for E in Energies ] )
            
            Position_Spacing = Parameters[ "Position Spacing (nm)" ]
                        
            Incoherent_Substrate = Parameters[ "Incoherent Substrate" ]
            Infinite_Final_Layer = Parameters[ "Infinite Final Layer" ] 
            
            Electric_Field_Output, Optical_Modelling_Output, Generation_Rate_Output = Single_Pass_Optical_Modeller( 
                
                Energies, Wavelengths, IQEs, Incident_Photon_Fluxes, 
                                             
                Surrounding_Medium, Layer_Types, Layer_Thicknesses, Position_Spacing,
                                 
                Active_Layer_Index, Incoherent_Substrate,Infinite_Final_Layer )
            
            EQEs = Optical_Modelling_Output[ "EQE" ]

            if Parameters[ "Junction Type" ] == "Multi-Junction":
                Bottom_Active_Layer_Index = Parameters[ "Bottom Active Layer Index" ] 
                BC_EQEs = Optical_Modelling_Output[ "Abs_" + str( Bottom_Active_Layer_Index ) + '_' + Layer_Types[ Bottom_Active_Layer_Index - 1 ] ]

            DataFrame( Optical_Modelling_Output ).to_excel( path.join( Output_Path, "Spectral_Output.xlsx" ), index = False )
            DataFrame( Electric_Field_Output ).to_excel( path.join( Output_Path, "Optical_E_Field.xlsx" ), index = False )
            DataFrame( Generation_Rate_Output ).to_excel( path.join( Output_Path, "Generation_Rate.xlsx" ), index = False )
                
        if Generate_Plots:
                
            plt.figure( dpi = 200 )
            if Parameters[ "Optically Modelled" ] == "True": 
                plt.plot( Energies, EQEs, label = 'EQE' )
                if Parameters[ "Junction Type" ] == "Multi-Junction":
                    plt.plot( Energies, BC_EQEs, label = 'Bottom Cell EQE' )
                plt.plot( Energies, Optical_Modelling_Output[ "R" ], label = 'R' )
                plt.plot( Energies, Optical_Modelling_Output[ "T" ], label = 'T' )
                plt.legend()
            else:
                plt.plot( Energies, EQEs )
            plt.xlabel( 'Energy, $E$ (eV)' )
            plt.ylabel( '$\\mathrm{EQE}_{\\mathrm{PV}}$' )
            plt.savefig( path.join( Output_Path, 'EQE.jpg' ) )      
    
            plt.figure( dpi = 200 )
            if Parameters[ "Optically Modelled" ] == "True": 
                plt.plot( Energies, EQEs, label = 'EQE' )
                if Parameters[ "Junction Type" ] == "Multi-Junction":
                    plt.plot( Energies, BC_EQEs, label = 'Bottom Cell EQE' )
                plt.plot( Energies, Optical_Modelling_Output[ "R" ], label = 'R' )
                plt.plot( Energies, Optical_Modelling_Output[ "T" ], label = 'T' )
                plt.legend()
            else:
                plt.plot( Energies, EQEs )
            plt.xlabel( 'Energy, $E$ (eV)' )
            plt.ylabel( '$\\mathrm{EQE}_{\\mathrm{PV}}$' )
            plt.yscale( 'log' )
            plt.savefig( path.join( Output_Path, 'EQE_logscale.jpg' ) )                  

##            if Parameters[ "Optically Modelled" ] == "True": 
                # Plot additional graphs
                
#                DataFrame( Optical_Modelling_Output ).to_excel( path.join( Output_Path, "Spectral_Output.xlsx" ), index = False )

    if Mode == "Many-System Comparison":
        
        # Experimental By Default; Use Cross-Referenced Filenames
               
        E_gs = { System : Systems_Database[ "E_g [eV]" ][ System ] for System in System_Identifiers }
        E_lowers = { System : Systems_Database[ "E_lower [eV]" ][ System ] for System in System_Identifiers }
        V_oc_exps = { System : Systems_Database[ 'V_oc [V]' ][ System ] for System in System_Identifiers }
    
        EQEs = { System : EQE_Loader_and_Interpolator( System + '.xlsx', Energies, E_lowers[ System ] ) 
                for System in System_Identifiers }
    
    if Mode == "Varied Lower Limit":
        
        EQE_Filename = Parameters[ "Experimental EQE Model" ]
        
        if '.xlsx' in EQE_Filename:
            
            EQE_Filename.replace( '.xlsx', '')
            
        E_g = Systems_Database[ "E_g [eV]" ][ EQE_Filename ]
        
        V_oc_exp = Systems_Database[ 'V_oc [V]' ][ EQE_Filename ]
        
        Data = read_excel( path.join( EQE_Spectra_Folder_Path, EQE_Filename + '.xlsx' ) )
        
        Wavelength_Lowers = array( Data[ Data.columns[ 0 ] ] ) 
        
        E_Lowers = Energy_Wavelength_Converter( Wavelength_Lowers )
        
        EQEs_initial = EQE_Loader_and_Interpolator( EQE_Filename + '.xlsx', Energies, 0 )
    
        if Generate_Plots:
            plt.figure( dpi = 200 )
            plt.plot( Energies, EQEs_initial, label = 'Top Cell' )
            plt.xlabel( '$E$ (eV)' )
            plt.ylabel( '$\\mathrm{EQE}_{\\mathrm{PV}}$' )
            plt.savefig( path.join( Output_Path, 'EQE.jpg' ) )
        
            plt.figure( dpi = 200 )
            plt.plot( Energies, EQEs_initial, label = 'Top Cell' )
            plt.xlabel( '$E$ (eV)' )
            plt.ylabel( '$\\mathrm{EQE}_{\\mathrm{PV}}$' )
            plt.yscale( 'log' )
            plt.savefig( path.join( Output_Path, 'EQE_logscale.jpg' ) )
            
            
        EQEs = { E_Lower : array( [ EQEs_initial[ j ] if E >= E_Lower else 0 for j, E in enumerate( Energies ) ] )
                
                for E_Lower in E_Lowers }
        
    if Mode == "Varied Gap":
        
        # Only done in a simualted simulation
        
        Simulated_EQE_Type = Parameters[ "Simulated EQE Model" ]
    
        Minimum_E_g = Parameters[ "Minimum E_g (eV)" ]
        Maximum_E_g = Parameters[ "Maximum E_g (eV)" ]
        N_E_g = Parameters[ "Number of E_gs" ]
    
        E_gs = linspace( Minimum_E_g, Maximum_E_g, N_E_g )
        
        EQE_max = Parameters[ "Above Gap EQE" ]
        
        if Simulated_EQE_Type == "Step Function":
        
            EQE_min = Parameters[ "Below Gap EQE" ]
    
            EQEs = { E_g : Step_Function_EQE_Simulator( Energies, E_g, EQE_min, EQE_max )
                    
                    for E_g in E_gs }
    
        if Simulated_EQE_Type == "Urbach Tail":
            
            E_U = Parameters[ "Urbach Energy (eV)" ]
    
            EQEs = { E_g : Urbach_Tail_EQE_Simulator( Energies, E_g, E_U, EQE_max )
                    
                    for E_g in E_gs }
    
        if Simulated_EQE_Type == "Disordered Organic PV":
            
            sigma_s = Parameters[ "Excitonic Disorder (eV)" ]
    
            EQEs = { E_g : Exction_Absorption_EQE_Simulator( Energies, E_g, sigma_s, EQE_max, kT )
                    
                    for E_g in E_gs }
                  
    # Multi-Junction with one Fixed Cell - Uses Elements from Above
    
    if Mode in [ "Multi-Junction - One Fixed Cell" ]:
        
        EQE_Mode = Parameters[ "Material Model" ] # Top Cell EQE Mode
        Bottom_Cell_EQE_Mode = Parameters[ "Bottom Cell Material Model" ]
        
        #### IF TOP CELL IS FIXED
        
        if Parameters[ "Fixed Multi-Junction Cell" ] == "Top" :
            
            # Top Cell is Fixed, Bottom is Varied; Load Top Cell EQEs
            
            #--------------------------------------------------------------------------
            if EQE_Mode == "Experimental (Loaded)":
            #--------------------------------------------------------------------------
                
                EQE_Filename = Parameters[ "Experimental EQE Model" ]
                
                if '.xlsx' in EQE_Filename:
                    
                    EQE_Filename.replace( '.xlsx', '')
                    
                E_g = Systems_Database[ "E_g [eV]" ][ EQE_Filename ]
                E_lower = Systems_Database[ "E_lower [eV]" ][ EQE_Filename ]
                V_oc_exp = Systems_Database[ 'V_oc [V]' ][ EQE_Filename ]
        
                EQEs = EQE_Loader_and_Interpolator( EQE_Filename + '.xlsx', Energies, E_lower )
            
                # print( EQE_Filename, ': E_g = ', E_g, 'eV; V_oc = ', V_oc_exp, 'V' )            
            
            #--------------------------------------------------------------------------
            if EQE_Mode == "Simulated":
            #--------------------------------------------------------------------------
                    
                Simulated_EQE_Type = Parameters[ "Simulated EQE Model" ]
            
                E_g = Parameters[ "E_g (eV)" ]
            
                EQE_max = Parameters[ "Above Gap EQE" ]
                
                if Simulated_EQE_Type == "Step Function":
                
                    EQE_min = Parameters[ "Below Gap EQE" ]
            
                    EQEs = Step_Function_EQE_Simulator( Energies, E_g, EQE_min, EQE_max )
            
                if Simulated_EQE_Type == "Urbach Tail":
                    
                    E_U = Parameters[ "Urbach Energy (eV)" ]
        
                    EQEs = Urbach_Tail_EQE_Simulator( Energies, E_g, E_U, EQE_max )
            
                if Simulated_EQE_Type == "Disordered Organic PV":
                    
                    sigma_s = Parameters[ "Excitonic Disorder (eV)" ]
        
                    EQEs = Exction_Absorption_EQE_Simulator( Energies, E_g, sigma_s, EQE_max, kT )
                
            # Load Bottom Cell EQEs
            
            #------------------------------------------------------------------
            if Bottom_Cell_EQE_Mode == "Experimental (Loaded)":
            #------------------------------------------------------------------
    
                BC_E_gs = { System : Systems_Database[ "E_g [eV]" ][ System ] for System in System_Identifiers }
                BC_E_lowers = { System : Systems_Database[ "E_lower [eV]" ][ System ] for System in System_Identifiers }
                BC_V_oc_exps = { System : Systems_Database[ 'V_oc [V]' ][ System ] for System in System_Identifiers }
    
                BC_EQEs = { System : EQE_Loader_and_Interpolator( System + '.xlsx', Energies, BC_E_lowers[ System ] ) 
                           for System in System_Identifiers }

            #------------------------------------------------------------------
            if Bottom_Cell_EQE_Mode == "Simulated":
            #------------------------------------------------------------------
                    
                BC_Simulated_EQE_Type = Parameters[ "Bottom Cell Simulated EQE Model" ]
            
                BC_E_gs = linspace( Parameters[ "Minimum E_g (eV)" ],
                                   
                                   Parameters[ "Maximum E_g (eV)" ],
                                   
                                   Parameters[ "Number of E_gs" ] )
                
                BC_EQEs = {}
            
                BC_EQE_max = Parameters[ "Bottom Cell Above Gap EQE" ]
                
                for BC_E_g in BC_E_gs:
                    
                    if BC_Simulated_EQE_Type == "Step Function":
                    
                        BC_EQE_min = Parameters[ "Bottom Cell Below Gap EQE" ]
                
                        BC_EQEs[ BC_E_g ] = Step_Function_EQE_Simulator( Energies, BC_E_g, BC_EQE_min, BC_EQE_max )
                
                    if BC_Simulated_EQE_Type == "Urbach Tail":
                        
                        BC_E_U = Parameters[ "Bottom Cell Urbach Energy (eV)" ]
            
                        BC_EQEs[ BC_E_g ] = Urbach_Tail_EQE_Simulator( Energies, BC_E_g, BC_E_U, BC_EQE_max )
                
                    if BC_Simulated_EQE_Type == "Disordered Organic PV":
                        
                        BC_sigma_s = Parameters[ "Bottom Cell Excitonic Disorder (eV)" ]
            
                        BC_EQEs[ BC_E_g ] = Exction_Absorption_EQE_Simulator( Energies, BC_E_g, BC_sigma_s, BC_EQE_max, kT )

        if Parameters[ "Fixed Multi-Junction Cell" ] == "Bottom" :
           
           # Top Cell is Fixed, Bottom is Varied; Load Bottom Cell EQEs
           
           #--------------------------------------------------------------------------
           if Bottom_Cell_EQE_Mode == "Experimental (Loaded)":
           #--------------------------------------------------------------------------
               
               BC_EQE_Filename = Parameters[ "Bottom Cell Experimental EQE Model" ]
               
               if '.xlsx' in BC_EQE_Filename:
                   
                   BC_EQE_Filename.replace( '.xlsx', '')
                   
               BC_E_g = Systems_Database[ "E_g [eV]" ][ BC_EQE_Filename ]
               BC_E_lower = Systems_Database[ "E_lower [eV]" ][ BC_EQE_Filename ]
               BC_V_oc_exp = Systems_Database[ 'V_oc [V]' ][ BC_EQE_Filename ]
       
               BC_EQEs = EQE_Loader_and_Interpolator( BC_EQE_Filename + '.xlsx', Energies, BC_E_lower )
                      
           #--------------------------------------------------------------------------
           if Bottom_Cell_EQE_Mode == "Simulated":
           #--------------------------------------------------------------------------
                   
               BC_Simulated_EQE_Type = Parameters[ "Simulated EQE Model" ]
           
               BC_E_g = Parameters[ "E_g (eV)" ]
           
               BC_EQE_max = Parameters[ "Above Gap EQE" ]
               
               if BC_Simulated_EQE_Type == "Step Function":
               
                   BC_EQE_min = Parameters[ "Below Gap EQE" ]
           
                   BC_EQEs = Step_Function_EQE_Simulator( Energies, BC_E_g, BC_EQE_min, BC_EQE_max )
           
               if BC_Simulated_EQE_Type == "Urbach Tail":
                   
                   BC_E_U = Parameters[ "Urbach Energy (eV)" ]
       
                   BC_EQEs = Urbach_Tail_EQE_Simulator( Energies, BC_E_g, BC_E_U, BC_EQE_max )
           
               if BC_Simulated_EQE_Type == "Disordered Organic PV":
                   
                   BC_sigma_s = Parameters[ "Excitonic Disorder (eV)" ]
       
                   BC_EQEs = Exction_Absorption_EQE_Simulator( Energies, BC_E_g, BC_sigma_s, BC_EQE_max, kT )
               
           # Load Top Cell EQEs
           
           #------------------------------------------------------------------
           if EQE_Mode == "Experimental (Loaded)":
           #------------------------------------------------------------------
   
               E_gs = { System : Systems_Database[ "E_g [eV]" ][ System ] for System in System_Identifiers }
               E_lowers = { System : Systems_Database[ "E_lower [eV]" ][ System ] for System in System_Identifiers }
               V_oc_exps = { System : Systems_Database[ 'V_oc [V]' ][ System ] for System in System_Identifiers }
   
               EQEs = { System : EQE_Loader_and_Interpolator( System + '.xlsx', Energies, E_lowers[ System ] ) 
                          for System in System_Identifiers }

           #------------------------------------------------------------------
           if EQE_Mode == "Simulated":
           #------------------------------------------------------------------
                   
               Simulated_EQE_Type = Parameters[ "Bottom Cell Simulated EQE Model" ]
           
               E_gs = linspace( Parameters[ "Minimum E_g (eV)" ],
                                  
                                  Parameters[ "Maximum E_g (eV)" ],
                                  
                                  Parameters[ "Number of E_gs" ] )
               
               EQEs = {}
           
               EQE_max  = Parameters[ "Bottom Cell Above Gap EQE" ]
               
               for E_g in E_gs:
                   
                   if Simulated_EQE_Type == "Step Function":
                   
                       EQE_min = Parameters[ "Bottom Cell Below Gap EQE" ]
               
                       EQEs[ E_g ] = Step_Function_EQE_Simulator( Energies, E_g, EQE_min, EQE_max )
               
                   if Simulated_EQE_Type == "Urbach Tail":
                       
                       E_U = Parameters[ "Bottom Cell Urbach Energy (eV)" ]
           
                       EQEs[ E_g ] = Urbach_Tail_EQE_Simulator( Energies, E_g, E_U, EQE_max  )
               
                   if Simulated_EQE_Type == "Disordered Organic PV":
                       
                       sigma_s = Parameters[ "Bottom Cell Excitonic Disorder (eV)" ]
           
                       EQEs[ E_g ] = Exction_Absorption_EQE_Simulator( Energies, E_g, sigma_s , EQE_max , kT )

    #==============================================================================
    # Part iii - Non-Radiative Open-Circuit Voltage Losses
    #==============================================================================
        
    if Mode in [ "Single Simulation", 'Varied Intensity', 'Varied Light Source', 'Typical Indoor Day' ]:
    
        # Use a single NR value in these cases (either radiative limit, empirically modelled, or determined from one-Sun Voc)
    
        NR_Loss_Mode = Parameters[ "Non-Radiative Loss Mode" ]
        
        if NR_Loss_Mode == "Radiative":
            EQE_EL = 1 
            Delta_V_oc_nr = 0
    
        if NR_Loss_Mode == "Empirical Organic PV Model":
            Delta_V_oc_nr = Delta_V_oc_NR_Empirical( E_g ) 
            EQE_EL = exp( - Delta_V_oc_nr / kT )        
                                    
        if NR_Loss_Mode == "Estimated from One-Sun Voc":
            
            EQE_EL_rad = 1
            
            Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                    
                               EQEs, 
                               
                               One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_rad,
                               
                               k * 298.15 / e,   # USE 25 deg C
                               
                               R_sh,
                               R_s,
                               
                               Area_cell )
            
            Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - V_oc_exp
            Delta_V_oc_nr = max( 0, Delta_V_oc_nr ) # ensure minimum of 0
            EQE_EL = exp( - Delta_V_oc_nr / kT )       
        
        if Multijunction_Simulation:
            
            #------------------------------------------------------------------
            # Assume Above NR Loss is for Top Cell; Det. Nr Loss for Bottom Cell
            #------------------------------------------------------------------ 
            
            BC_NR_Loss_Mode = Parameters[ "Bottom Cell Non-Radiative Loss Mode" ]
            
            if BC_NR_Loss_Mode == "Radiative":
                BC_EQE_EL = 1 
                BC_Delta_V_oc_nr = 0
        
            if BC_NR_Loss_Mode == "Empirical Organic PV Model":
                BC_Delta_V_oc_nr = Delta_V_oc_NR_Empirical( BC_E_g ) 
                BC_EQE_EL = exp( - BC_Delta_V_oc_nr / kT )        
                                        
            if BC_NR_Loss_Mode == "Estimated from One-Sun Voc":
                EQE_EL_rad = 1
                
                BC_Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                        
                                   BC_EQEs, 
                                   
                                   One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_rad,
                                   
                                   k * 298.15 / e,   # USE 25 deg C
                                   
                                   BC_R_sh,
                                   BC_R_s,
                                   
                                   Area_cell )
                
                BC_Delta_V_oc_nr = BC_Radiative_Output[ 'V_oc' ] - BC_V_oc_exp
                BC_Delta_V_oc_nr = max( 0, BC_Delta_V_oc_nr ) # ensure minimum of 0
                BC_EQE_EL = exp( - BC_Delta_V_oc_nr / kT )       

    if Mode == "Many-System Comparison":
    
        EQE_EL_rad = 1
    
        Delta_V_oc_nrs = {}
        EQE_ELs = {}
        
        for System in System_Identifiers:
    
            Radiative_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                               EQEs[ System ], 
                               
                               One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_rad,
                               
                               k * 298.15 / e, # USE 25 deg C
                               
                               R_sh,
                               R_s,
                               
                               Area_cell )
    
            Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - V_oc_exps[ System ] 
    
            Delta_V_oc_nr = max( 0, Delta_V_oc_nr ) # ensure minimum of 0
            
            EQE_EL = exp( - Delta_V_oc_nr / kT )
    
            Delta_V_oc_nrs[ System ] = Delta_V_oc_nr
            EQE_ELs[ System ] = EQE_EL        
        
    if Mode == "Varied Lower Limit":
        
        V_oc_rads = {}
        Delta_V_oc_nrs = {}
        EQE_ELs = {}
        
        EQE_EL_rad = 1
        
        for E_Lower in E_Lowers: 
            
            Radiative_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                               EQEs[ E_Lower ], 
                               
                               One_Sun_Photon_Fluxes,     #### USE ONE-SUN FLUXES
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL_rad,
                               
                               k * 298.15 / e,  # USE 25 deg C
                               
                               R_sh,
                               R_s,
                               
                               Area_cell )
            
            V_oc_rads[ E_Lower ] = Radiative_Output[ 'V_oc' ]
            Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - V_oc_exp
            Delta_V_oc_nr = max( 0, Delta_V_oc_nr ) # ensure minimum of 0
            EQE_EL = exp( - Delta_V_oc_nr / kT )        
        
            Delta_V_oc_nrs[ E_Lower ] = Delta_V_oc_nr
            EQE_ELs[ E_Lower ] = EQE_EL
        
    if Mode == "Varied Gap":
        
        NR_Loss_Mode = Parameters[ "Non-Radiative Loss Mode" ]
    
        if NR_Loss_Mode == 'Radiative':
            
            Delta_V_oc_nrs = { E_g : 0 for E_g in E_gs }
            EQE_ELs = { E_g : 1 for E_g in E_gs }
            
        if NR_Loss_Mode == "Empirical Organic PV Model":
            
            Delta_V_oc_nrs = { E_g : Delta_V_oc_NR_Empirical( E_g ) for E_g in E_gs }
            EQE_ELs = { E_g : exp( - Delta_V_oc_nrs[ E_g ] / kT ) for E_g in E_gs }
            
        if NR_Loss_Mode == 'Estimated from One-Sun Voc':
            
            print( "ERROR - Update NR Loss model for simulated EQE. Use either (i) Radiative or (ii) Empirical Organic PV Model." )
    
            exit()
            
    if Mode == "Multi-Junction - One Fixed Cell":
        
        #----------------------------------------------------------------------        
        if Parameters[ "Fixed Multi-Junction Cell" ] == "Top":
            
            # Determine Top Cell NR Loss then Bottom Cell NR Losses
        
            NR_Loss_Mode = Parameters[ "Non-Radiative Loss Mode" ]
            
            if NR_Loss_Mode == "Radiative":
                EQE_EL = 1 
                Delta_V_oc_nr = 0
        
            if NR_Loss_Mode == "Empirical Organic PV Model":
                Delta_V_oc_nr = Delta_V_oc_NR_Empirical( E_g ) 
                EQE_EL = exp( - Delta_V_oc_nr / kT )        
                                        
            if NR_Loss_Mode == "Estimated from One-Sun Voc":
                
                EQE_EL_rad = 1
                
                Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                        
                                   EQEs, 
                                   
                                   One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   EQE_EL_rad,
                                   
                                   k * 298.15 / e,   # USE 25 deg C
                                   
                                   R_sh,
                                   R_s,
                                   
                                   Area_cell )
                
                Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - V_oc_exp
                Delta_V_oc_nr = max( 0, Delta_V_oc_nr ) # ensure minimum of 0
                EQE_EL = exp( - Delta_V_oc_nr / kT ) 
            #------------------------------------------------------------------
            # Evaluate Bottom Cell NR Losses
            
            BC_NR_Loss_Mode = Parameters[ "Non-Radiative Loss Mode" ]
            
            if BC_NR_Loss_Mode == "Radiative":
                BC_EQE_ELs = { Key : 1 for Key in list( BC_EQEs.keys() ) } 
                BC_Delta_V_oc_nrs = { Key : 0 for Key in list( BC_EQEs.keys() ) } 
        
            if BC_NR_Loss_Mode == "Empirical Organic PV Model":
                BC_Delta_V_oc_nrs = { Key : Delta_V_oc_NR_Empirical( BC_E_gs[ Key] ) for Key in list( BC_EQEs.keys() ) } 
                BC_EQE_ELs = { Key : exp( - BC_Delta_V_oc_nrs[ Key ] / kT ) for Key in list( BC_EQEs.keys() ) }
                                        
            if BC_NR_Loss_Mode == "Estimated from One-Sun Voc":
                
                BC_EQE_EL_rad = 1
                BC_Delta_V_oc_nrs = {}
                BC_EQE_ELs = {}
                
                for Key in list( BC_EQEs.keys() ): 
                
                    Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                            
                                       BC_EQEs[ Key ], 
                                       
                                       One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                                       
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
                                       
                                       BC_EQE_EL_rad,
                                       
                                       k * 298.15 / e,   # USE 25 deg C
                                       
                                       BC_R_sh,
                                       BC_R_s,
                                       
                                       Area_cell )
                    
                    BC_Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - BC_V_oc_exps[ Key ]
                    BC_Delta_V_oc_nr = max( 0, BC_Delta_V_oc_nr ) # ensure minimum of 0
                    BC_EQE_EL = exp( - BC_Delta_V_oc_nr / kT )    
            
                    BC_Delta_V_oc_nrs[ Key ] = BC_Delta_V_oc_nr
                    BC_EQE_ELs[ Key ] = BC_EQE_EL
                    
        #----------------------------------------------------------------------            
        if Parameters[ "Fixed Multi-Junction Cell" ] == "Bottom":
            
            # Determine Bottom Cell NR Loss then Top Cell NR Losses
        
            BC_NR_Loss_Mode = Parameters[ "Bottom Cell Non-Radiative Loss Mode" ]
            
            if BC_NR_Loss_Mode == "Radiative":
                BC_EQE_EL = 1 
                BC_Delta_V_oc_nr = 0
        
            if BC_NR_Loss_Mode == "Empirical Organic PV Model":
                BC_Delta_V_oc_nr = Delta_V_oc_NR_Empirical( E_g ) 
                BC_EQE_EL = exp( - Delta_V_oc_nr / kT )        
                                        
            if BC_NR_Loss_Mode == "Estimated from One-Sun Voc":
                
                BC_EQE_EL_rad = 1
                
                Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                        
                                   BC_EQEs, 
                                   
                                   One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,
                                   
                                   BC_EQE_EL_rad,
                                   
                                   k * 298.15 / e,   # USE 25 deg C
                                   
                                   BC_R_sh,
                                   BC_R_s,
                                   
                                   Area_cell )
                
                BC_Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - BC_V_oc_exp
                BC_Delta_V_oc_nr = max( 0, BC_Delta_V_oc_nr ) # ensure minimum of 0
                BC_EQE_EL = exp( - BC_Delta_V_oc_nr / kT )    
                
            #------------------------------------------------------------------
            # Evaluate Top Cell NR Losses            
            
            NR_Loss_Mode = Parameters[ "Non-Radiative Loss Mode" ]
            
            if NR_Loss_Mode == "Radiative":
                EQE_ELs = { Key : 1 for Key in list( EQEs.keys() ) } 
                Delta_V_oc_nrs = { Key : 0 for Key in list( EQEs.keys() ) } 
        
            if NR_Loss_Mode == "Empirical Organic PV Model":
                Delta_V_oc_nrs  = { Key : Delta_V_oc_NR_Empirical( BC_E_gs[ Key] ) for Key in list( EQEs.keys() ) } 
                EQE_ELs = { Key : exp( - Delta_V_oc_nrs[ Key ] / kT ) for Key in list( EQEs.keys() ) }
                                        
            if NR_Loss_Mode == "Estimated from One-Sun Voc":
                
                EQE_EL_rad = 1
                Delta_V_oc_nrs  = {}
                EQE_ELs = {}
                
                for Key in list( EQEs.keys() ): 
                
                    Radiative_Output = Small_Area_Figure_of_Merit_Calculator( 100.3, # one-Sun Intensity
                            
                                       EQEs[ Key ], 
                                       
                                       One_Sun_Photon_Fluxes,  # USE ONE-SUN FLUXES
                                       
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
                                       
                                       EQE_EL_rad,
                                       
                                       k * 298.15 / e,   # USE 25 deg C
                                       
                                       BC_R_sh,
                                       BC_R_s,
                                       
                                       Area_cell )
                    
                    Delta_V_oc_nr = Radiative_Output[ 'V_oc' ] - V_oc_exps[ Key ]
                    Delta_V_oc_nr  = max( 0, Delta_V_oc_nr  ) # ensure minimum of 0
                    EQE_EL = exp( - Delta_V_oc_nr  / kT )    
            
                    Delta_V_oc_nrs [ Key ] = Delta_V_oc_nr 
                    EQE_ELs[ Key ] = EQE_EL
                    
    #==========================================================================
    # Part iv - Run Small-Area Simulations 
    #==========================================================================
        
    global Output
    
    Output = {}
    
    if Multijunction_Simulation:
        
        Terminal_Type = Parameters[ "Multi-Junction Terminal Type" ]
        Coupling_1_to_2 = Parameters[ "Radiative Coupling Top to Bottom" ]
        Coupling_2_to_1 = Parameters[ "Radiative Coupling Bottom to Top" ]
        
        if Coupling_1_to_2 == 'True':
            Coupling_1_to_2 = True
        else:
            Coupling_1_to_2 = False
        
        if Coupling_2_to_1 == 'True':
            Coupling_2_to_1 = True
        else:
            Coupling_2_to_1 = False
        
        
    #--------------------------------------------------------------------------
    # 1. Single Simulation
    #--------------------------------------------------------------------------
            
    if Mode == "Single Simulation":

        if not Multijunction_Simulation:
            
            #------------------------------------------------------------------
            # Treat Single Junction Case First
            #------------------------------------------------------------------
            
            Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                               EQEs, 
                               
                               Incident_Photon_Fluxes,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_EL,
                               
                               kT,
                               
                               R_sh,
                               R_s,
                               
                               Area_cell )
            
            Single_Output[ 'EQE_EL' ] = EQE_EL
            Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nr 
            
        if Multijunction_Simulation:

            #------------------------------------------------------------------
            # Treat Multi-Junction Case First
            #------------------------------------------------------------------
                        
            Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                               EQEs, 
                               BC_EQEs, 
                               
                               Incident_Photon_Fluxes,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,

                               EQE_EL,
                               BC_EQE_EL,                         
                               
                               kT,
                               
                               R_sh,
                               BC_R_sh,
                                            
                               R_s,
                               BC_R_s,
                               Area_cell,

                               Terminal_Type = Terminal_Type,
                               
                               Coupling_1_to_2 = Coupling_1_to_2,
                               Coupling_2_to_1 = Coupling_2_to_1,
                               
                               Convergence_Criterion = 1e-6, 
                               Max_Iterations = 300 )
            
            Single_Output[ 'Top Cell EQE_EL' ] = EQE_EL
            Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_EL
            
            Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nr 
            Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nr 
        
        Output = { Key : [ Single_Output[ Key ] ] for Key in list( Single_Output.keys() ) }
    
    #--------------------------------------------------------------------------
    # 2. Many-System Comparison
    #--------------------------------------------------------------------------
                
    if Mode == "Many-System Comparison":
        
        Output[ 'Systems' ] = System_Identifiers
        
        Output[ 'E_g' ] = list( E_gs.values() )
        
        Output[ 'E_lower' ] = list( E_lowers.values() )
        
        Output[ 'V_oc_exp' ] = list( V_oc_exps.values() )
        
        for j, System in enumerate( System_Identifiers ):
    
            Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                               EQEs[ System ], 
                               
                               Incident_Photon_Fluxes,
                               
                               Black_Body_Fluxes,
                               
                               Energies,
                               
                               Wavelengths,
                               
                               EQE_ELs[ System ],
                               
                               kT,
                               
                               R_sh,
                               R_s,
                               
                               Area_cell )
        
            Single_Output[ 'EQE_EL' ] = EQE_ELs[ System ]
        
            Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nrs[ System ] 
        
            if j == 0:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ] = [ Single_Output[ Key ] ]
                    
            else:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ].append( Single_Output[ Key ] )
        
    #--------------------------------------------------------------------------
    # 3. Varied Gap
    #--------------------------------------------------------------------------
              
    if Mode == "Varied Gap":
    
        Output[ 'E_g' ] = list( E_gs )
        
        for j, E_g in enumerate( E_gs ):
    
            Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                                EQEs[ E_g ], 
                               
                                Incident_Photon_Fluxes,
                               
                                Black_Body_Fluxes,
                               
                                Energies,
                               
                                Wavelengths,
                               
                                EQE_ELs[ E_g ],
                               
                                kT,
                               
                                R_sh,
                                R_s,
                               
                                Area_cell )
        
            Single_Output[ 'EQE_EL' ] = EQE_ELs[ E_g ]
        
            Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nrs[ E_g ]
        
            if j == 0:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ] = [ Single_Output[ Key ] ]
                    
            else:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ].append( Single_Output[ Key ] )
        
    #--------------------------------------------------------------------------
    # 4. Varied Lower Limit
    #--------------------------------------------------------------------------
                    
    if Mode == "Varied Lower Limit":
        
        Output[ 'E_Lower' ] = list( E_Lowers )
        
        for j, E_Lower in enumerate( E_Lowers ):
    
            Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                                EQEs[ E_Lower ], 
                               
                                Incident_Photon_Fluxes,
                               
                                Black_Body_Fluxes,
                               
                                Energies,
                               
                                Wavelengths,
                               
                                EQE_ELs[ E_Lower ],
                               
                                kT,
                               
                                R_sh,
                                R_s,
                               
                                Area_cell )
        
            Single_Output[ 'V_oc_rad' ] = V_oc_rads[ E_Lower ]
            
            Single_Output[ 'EQE_EL' ] = EQE_ELs[ E_Lower ]
        
            Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nrs[ E_Lower ]
        
            if j == 0:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ] = [ Single_Output[ Key ] ]
                    
            else:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ].append( Single_Output[ Key ] )
        
    #--------------------------------------------------------------------------
    # 5. Varied Intensity
    #--------------------------------------------------------------------------
             
    if Mode == "Varied Intensity":
        
        Output[ 'Intensity [mW/cm2]'] = list( Light_Powers )
        
        if Intensity_Unit != 'mw/cm2':
        
            Output[ 'Intensity [' + Intensity_Unit + ']' ] = list( Light_Powers_initial )
            
        for j, Light_Power in enumerate( Light_Powers ):
            
            print( Light_Power )
            
            if not Multijunction_Simulation:

                Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Power,
                    
                                EQEs, 
                               
                                Incident_Photon_Fluxes[ Light_Power ],
                               
                                Black_Body_Fluxes,
                               
                                Energies,
                               
                                Wavelengths,
                               
                                EQE_EL,
                               
                                kT,
                               
                                R_sh,
                                R_s,
                               
                                Area_cell )
        
                Single_Output[ 'EQE_EL' ] = EQE_EL
        
                Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nr
                
            if Multijunction_Simulation:
                
                Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Power,
                        
                                   EQEs, 
                                   BC_EQEs, 
                                   
                                   Incident_Photon_Fluxes[ Light_Power ],
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,

                                   EQE_EL,
                                   BC_EQE_EL,                         
                                   
                                   kT,
                                   
                                   R_sh,
                                   BC_R_sh,       
                                   R_s,
                                   BC_R_s,
                                   Area_cell,

                                   Terminal_Type = Terminal_Type,
                                   
                                   Coupling_1_to_2 = Coupling_1_to_2,
                                   Coupling_2_to_1 = Coupling_2_to_1,
                                   
                                   Convergence_Criterion = 1e-6, 
                                   Max_Iterations = 300 )
                
                Single_Output[ 'Top Cell EQE_EL' ] = EQE_EL
                Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_EL
                
                Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nr 
                Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nr 
        
            if j == 0:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ] = [ Single_Output[ Key ] ]
                    
            else:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ].append( Single_Output[ Key ] )
              
    #--------------------------------------------------------------------------
    # 6. Varied Light Source
    #--------------------------------------------------------------------------

    if Mode == "Varied Light Source":
        
        Output[ 'Light Source' ] = Select_Light_Sources
        
        Output[ 'Intensity [mW/cm2]' ] = list( Light_Powers.values() )
                
        if Intensity_Unit != 'mw/cm2':
        
            Output[ 'Intensity [' + Intensity_Unit + ']' ] = [ Light_Powers_initial[ Light_Source ] for Light_Source in Select_Light_Sources ]
            
        for j, Light_Source in enumerate( Select_Light_Sources ):
            
            I_j = Light_Powers[ Light_Source ] 
            
            print( Light_Source, I_j )
            
            if I_j != 0:
                
                # Only Simulate performance if irradiance is not zero
                
                if not Multijunction_Simulation:
            
                    Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Powers[ Light_Source ],
                            
                                       EQEs, 
                                       
                                       Incident_Photon_Fluxes[ Light_Source ],
                                       
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
                                       
                                       EQE_EL,
                                       
                                       kT,
                                       
                                       R_sh,
                                       R_s,
                                       
                                       Area_cell )
                
                    Single_Output[ 'EQE_EL' ] = EQE_EL
                
                    Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nr
                    
                if Multijunction_Simulation:
                    
                    Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Powers[ Light_Source ],
                            
                                       EQEs, 
                                       BC_EQEs, 
                                       
                                       Incident_Photon_Fluxes[ Light_Source ],
                                       
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
    
                                       EQE_EL,
                                       BC_EQE_EL,                         
                                       
                                       kT,
                                       
                                       R_sh,
                                       BC_R_sh,
                                                    
                                       R_s,
                                       BC_R_s,
                                       Area_cell,
    
                                       Terminal_Type = Terminal_Type,
                                       
                                       Coupling_1_to_2 = Coupling_1_to_2,
                                       Coupling_2_to_1 = Coupling_2_to_1,
                                       
                                       Convergence_Criterion = 1e-6, 
                                       Max_Iterations = 300 )
                    
                    Single_Output[ 'Top Cell EQE_EL' ] = EQE_EL
                    Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_EL
                    
                    Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nr 
                    Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nr 
                
                if j == 0:
                    
                    for Key in list( Single_Output.keys() ):
                                        
                        Output[ Key ] = [ Single_Output[ Key ] ]
                        
                else:
                    
                    for Key in list( Single_Output.keys() ):
                                        
                        Output[ Key ].append( Single_Output[ Key ] )
        
    #--------------------------------------------------------------------------
    # 7. Typical Indoor Day
    #--------------------------------------------------------------------------
        
    if Mode == 'Typical Indoor Day':    
        
        Times = [ int( Time[ 5: ] ) / 100 for Time in Day_Light_Sources ]
        
        Output[ "Time" ] = Times
        Output[ "Light Source" ] = Day_Light_Sources
        Output[ 'Intensity [mW/cm2]'] = list( Light_Powers.values() )
        
        if Intensity_Unit == 'W/m2':
            Output[ 'Intensity [W/m2]' ] = array( list( Light_Powers.values() ) ) * 10
        
        if Intensity_Unit == 'lux':
            Output[ 'Intensity [lux]' ] = [ 10 * Light_Powers[ Light_Source ] * Lux_Constant_Calculator( Light_Source, Wavelengths, Energies )
                                           for Light_Source in Day_Light_Sources ]
        
        for j, Light_Source in enumerate( Day_Light_Sources ):
            
            I_j = Light_Powers[ Light_Source ] 
            
            print( Light_Source, I_j )
            
            if I_j != 0:
                
                # Only Simulate performance if irradiance is not zero
    
                if not Multijunction_Simulation:
                    
                    Single_Output = Small_Area_Figure_of_Merit_Calculator( Light_Powers[ Light_Source ],
                            
                                        EQEs, 
                                       
                                        Incident_Photon_Fluxes[ Light_Source ],
                                       
                                        Black_Body_Fluxes,
                                       
                                        Energies,
                                       
                                        Wavelengths,
                                       
                                        EQE_EL,
                                       
                                        kT,
                                       
                                        R_sh,
                                        R_s,
                                       
                                        Area_cell )
                    
                    Single_Output[ 'EQE_EL' ] = EQE_EL
                
                    Single_Output[ 'Delta_V_oc_nr' ] = Delta_V_oc_nr
            
                if Multijunction_Simulation:
                    
                    Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Powers[ Light_Source ],
                            
                                       EQEs, 
                                       BC_EQEs, 
                                       
                                       Incident_Photon_Fluxes[ Light_Source ],
                                       
                                       Black_Body_Fluxes,
                                       
                                       Energies,
                                       
                                       Wavelengths,
    
                                       EQE_EL,
                                       BC_EQE_EL,                         
                                       
                                       kT,
                                       
                                       R_sh,
                                       BC_R_sh,
                                                    
                                       R_s,
                                       BC_R_s,
                                       Area_cell,
    
                                       Terminal_Type = Terminal_Type,
                                       
                                       Coupling_1_to_2 = Coupling_1_to_2,
                                       Coupling_2_to_1 = Coupling_2_to_1,
                                       
                                       Convergence_Criterion = 1e-6, 
                                       Max_Iterations = 300 )
                    
                    Single_Output[ 'Top Cell EQE_EL' ] = EQE_EL
                    Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_EL
                    
                    Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nr 
                    Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nr 
        
            if j == 0:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ] = [ Single_Output[ Key ] ]
                    
            else:
                
                for Key in list( Single_Output.keys() ):
                                    
                    Output[ Key ].append( Single_Output[ Key ] )

    #--------------------------------------------------------------------------
    # 8. Multi-Junction - One Fixed Cell
    #--------------------------------------------------------------------------
                
    if Mode == "Multi-Junction - One Fixed Cell":
        
        if Parameters[ "Fixed Multi-Junction Cell" ] == "Top":
            
            if Parameters[ "Bottom Cell Material Model" ] == "Simulated":
                Output[ "Bottom Cell E_g (eV)" ] = list( BC_E_gs )
            else:
                Output[ "Bottom Cell E_g (eV)" ] = list( BC_E_gs.values() ) 
                Output[ "Bottom Cell System Type" ] = list( BC_E_gs.keys() ) 
                
            for j, Key in enumerate( list( BC_EQEs.keys() ) ):
                
                print( j, Key )
        
                Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Power,
                        
                                   EQEs, 
                                   BC_EQEs[ Key ], 
                                   
                                   Incident_Photon_Fluxes,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,

                                   EQE_EL,
                                   BC_EQE_ELs[ Key ],                         
                                   
                                   kT,
                                   
                                   R_sh,
                                   BC_R_sh,
                                                
                                   R_s,
                                   BC_R_s,
                                   Area_cell,

                                   Terminal_Type = Terminal_Type,
                                   
                                   Coupling_1_to_2 = Coupling_1_to_2,
                                   Coupling_2_to_1 = Coupling_2_to_1,
                                   
                                   Convergence_Criterion = 1e-6, 
                                   Max_Iterations = 300 )
                
                Single_Output[ 'Top Cell EQE_EL' ] = EQE_EL
                Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_ELs[ Key ]
                
                Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nr 
                Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nrs[ Key ]
    
                if j == 0:
                    
                    for Output_Key in list( Single_Output.keys() ):
                                        
                        Output[ Output_Key ] = [ Single_Output[ Output_Key ] ]
                        
                else:
                    
                    for Output_Key in list( Single_Output.keys() ):
                                        
                        Output[ Output_Key ].append( Single_Output[ Output_Key ] )
                    
        if Parameters[ "Fixed Multi-Junction Cell" ] == "Bottom":
                        
            if Parameters[ "Material Model" ] == "Simulated":
                
                Output[ "Top Cell E_g (eV)" ] = list( E_gs )
                
            
            else:
                
                Output[ "Top Cell E_g (eV)" ] = list( E_gs.values() ) 
                Output[ "Top Cell System Type" ] = list( E_gs.keys() ) 
                
            
            for j, Key in enumerate( list( EQEs.keys() ) ):
            
                Single_Output = MJ_Small_Area_Figure_of_Merit_Calculator( Light_Power,
                        
                                   EQEs[ Key ], 
                                   BC_EQEs, 
                                   
                                   Incident_Photon_Fluxes,
                                   
                                   Black_Body_Fluxes,
                                   
                                   Energies,
                                   
                                   Wavelengths,

                                   EQE_ELs[ Key ],
                                   BC_EQE_EL,                         
                                   
                                   kT,
                                   
                                   R_sh,
                                   BC_R_sh,
                                                
                                   R_s,
                                   BC_R_s,
                                   Area_cell,

                                   Terminal_Type = Terminal_Type,
                                   
                                   Coupling_1_to_2 = Coupling_1_to_2,
                                   Coupling_2_to_1 = Coupling_2_to_1,
                                   
                                   Convergence_Criterion = 1e-6, 
                                   Max_Iterations = 300 )
                
                Single_Output[ 'Top Cell EQE_EL' ] = EQE_ELs[ Key ]
                Single_Output[ 'Bottom Cell EQE_EL' ] = BC_EQE_EL
                
                Single_Output[ 'Top Cell Delta_V_oc_nr' ] = Delta_V_oc_nrs[ Key ]
                Single_Output[ 'Bottom Cell Delta_V_oc_nr' ] = BC_Delta_V_oc_nr
    
                if j == 0:
                    
                    for Output_Key in list( Single_Output.keys() ):
                                        
                        Output[ Output_Key ] = [ Single_Output[ Output_Key ] ]
                        
                else:
                    
                    for Output_Key in list( Single_Output.keys() ):
                                        
                        Output[ Output_Key ].append( Single_Output[ Output_Key ] )                    

    #==============================================================================
    # Part v - Run Large-Area Simulations (if True)
    #==============================================================================
    
    if Module_Simulation == 'True':
            
        Output_Copy = DataFrame( Output.copy() )
    
        #--------------------------------------------------------------------------
        # 1. Single Simulation
        #--------------------------------------------------------------------------
        
        if Mode == "Single Simulation":
            
            Large_Area_Output = Large_Area_FoM_Calculator( dict( Output_Copy.iloc[ 0 ] ),
                                                          
                           Light_Power,
                                      
                           Area_Module, Area_cell, N_Cells,
                           
                           L_x, L_y, delta,
                                              
                           Evaluate_System_Efficiency = True )
    
            Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                
            for Key in Additional_Keys:
                    
                Output[ Key ] = [ Large_Area_Output[ Key ] ]
                
        #--------------------------------------------------------------------------        
        # 2. Bulk Analysis            
        #--------------------------------------------------------------------------                    
    
        if Mode == "Many-System Comparison":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                            
                Large_Area_Output = Large_Area_FoM_Calculator( dict( Output_Copy.iloc[ Index ] ),
                                                              
                               Light_Power,
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                                                  
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                        
                    # print( Index )
                    
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )
            
        #--------------------------------------------------------------------------        
        # 3. Varied Light Source
        #--------------------------------------------------------------------------                    
    
        if Mode == "Varied Light Source":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Light_Powers[ Output_i[ "Light Source" ] ],
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )
                        
        #--------------------------------------------------------------------------        
        # 4. Varied Intensity
        #--------------------------------------------------------------------------                    
    
        if Mode == "Varied Intensity":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Light_Powers[ Index ],
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )                    
                        
        #--------------------------------------------------------------------------        
        # 5. Varied Gap
        #--------------------------------------------------------------------------                    
    
        if Mode == "Varied Gap":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Light_Power,
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )                    
                                            
        #--------------------------------------------------------------------------        
        # 6. Varied Lower Limit
        #--------------------------------------------------------------------------                    
    
        if Mode == "Varied Lower Limit":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Light_Power,
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )   
                                              
        #--------------------------------------------------------------------------        
        # 7. Typical Indoor Day
        #--------------------------------------------------------------------------                    
    
        if Mode == "Typical Indoor Day":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Output_i[ "Intensity [mW/cm2]" ],
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )         
                        
            # For the Typical Indoor Day Scenario, also determine energies.
            
            # Start by determining the incident energies:
                
            Incident_Energies = cumulative_simpson( Output[ "Intensity [mW/cm2]" ],
                                         
                                          x = Output[ "Time" ] ) # Hours, energy has units of mW hr/cm2
            
            Incident_Energies = Output[ "Area_Module" ] * array( [ 0 ] + list( Incident_Energies ) ) # add initial data point back
                        
            # Next, determine the converted energies:
            
            Converted_Energies = cumulative_simpson( Output[ "P_mpp_Module" ],
                                         
                                          x = Output[ "Time" ] ) # units of mW hr
            
            Converted_Energies = array( [ 0 ] + list( Converted_Energies ) ) # add initial data point back
    
            # Finally, determine the stored energies:
                
            Stored_Energies = cumulative_simpson( array( Output[ "P_mpp_Module" ] ) * array( Output[ "eta_interp" ],),
                                             
                                              x = Output[ "Time" ] ) # units of mW hr
                                
            Stored_Energies = array( [ 0 ] + list( Stored_Energies ) ) # add initial data point back
            
            Output = Output | { "Total Incident Energy (mWhr)" : Incident_Energies, 
                               
                               'Converted Energy (mWhr)' : Converted_Energies,
                               
                               'Stored Energy (mWhr)' : Stored_Energies }
            
        #--------------------------------------------------------------------------        
        # 8. Multi-Junction - One Fixed Cell
        #--------------------------------------------------------------------------                    
    
        if Mode == "One Fixed Cell":
    
            for Index in range( len( Output[ list( Output.keys() )[ 0 ] ] ) ):
                
                Output_i = dict( Output_Copy.iloc[ Index ] )
                            
                Large_Area_Output = Large_Area_FoM_Calculator( Output_i,
                                                              
                               Light_Power,
                                          
                               Area_Module, Area_cell, N_Cells,
                               
                               L_x, L_y, delta,
                   
                               Evaluate_System_Efficiency = True )
            
                Additional_Keys = [ Key for Key in list( Large_Area_Output.keys() ) if Key not in list( Output_Copy.keys() ) ]
                    
                for Key in Additional_Keys:
                                        
                    if Index == 0:
                                        
                        # Add additional keys
                        
                        Output[ Key ] = [ Large_Area_Output[ Key ] ]
                            
                    else:
                        
                        # print( Index )
                    
                        Output[ Key ].append( Large_Area_Output[ Key ] )               
            
    ### SAVE DATA
    
    DataFrame( Output ).to_csv( path.join( Output_Path, 'Output.csv' ), index = False )

    #==========================================================================
    # Generate Output Graphs
    #==========================================================================
    
    if Generate_Plots:
        
        # if Mode != "Single Simulation":
            
            # for all other modes
            
            # Independent_Variable = list( Output.keys() )[ 0 ]
            
            # Graphs = list( Output.keys() )[ 1: ]
            
            # for Graph in Graphs:
                
            #     plt.figure( dpi = 200 )
                
            #     plt.plot( Output[ Independent_Variable ], Output[ Graph ] )
                
            #     plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )

#------------------------------------------------------------------------------

        if Mode == "Many-System Comparison":
            
            E_gs = Output[ 'E_g' ]            
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'E_g' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.plot( E_gs, Output[ Graph ], '.' )
                    
                    plt.xlabel( Label_Dictionary[ 'E_g' ] )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                        
                    plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )

#------------------------------------------------------------------------------                    

        if Mode == "Varied Gap":
            
            E_gs = Output[ 'E_g' ]            
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'E_g' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.plot( E_gs, Output[ Graph ] )
                    
                    plt.xlabel( Label_Dictionary[ 'E_g' ] )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                        
                    plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )

#------------------------------------------------------------------------------
                    
        if Mode == "Varied Lower Limit":
            
            E_lowers = Output[ 'E_Lower' ]            
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'E_Lower' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.plot( E_lowers, Output[ Graph ] )
                    
                    plt.xlabel( Label_Dictionary[ 'E_lower' ] )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                        
                    plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )

#------------------------------------------------------------------------------

        if Mode == "Varied Intensity":
            
            Intensities = Output[ 'Intensity [mW/cm2]' ]            
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'Intensity [mW/cm2]' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.plot( Intensities, abs( array( Output[ Graph ] ) ) )
                    
                    plt.xlabel( Label_Dictionary[ 'Intensity [mW/cm2]' ] )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph in [ 'EQE_EL', 'J_sc', 'J_mpp', 'P_mpp', 'I_sc_Module', 'I_mpp_Module', 'P_mpp_Module' ]:
                        
                        plt.yscale( 'log' ) 
                        
                    if Parameters[ "Light Intensity Sample Type" ] == "Logarithmic":
                        
                        plt.xscale( 'log' )
                        
                    plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )
                    
#------------------------------------------------------------------------------

        if Mode == "Varied Light Source":
            
            Sources = Output[ 'Light Source' ]            
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'Light Source' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.bar( Sources, Output[ Graph ] )
                    
                    plt.xticks( rotation = 90 )
                    
                    plt.xlabel( 'Light Source' )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                    
                    if Graph != 'Intensity [mW/cm2]':

                        plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )
                    
#------------------------------------------------------------------------------

        if Mode == "Typical Indoor Day":
            
            Sources = Output[ "Time" ]       
                        
            Graphs = list( Output.keys() )
            
            Graphs.remove( 'Light Source' )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.plot( Sources, Output[ Graph ] )
                                        
                    plt.xlabel( 'Time of Day (Hours)' )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                    
                    if Graph != 'Intensity [mW/cm2]':

                        plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )
                    
#------------------------------------------------------------------------------

        if Mode == "Multi-Junction - One Fixed Cell":
            
            Graphs = list( Output.keys() )
            
            for Graph in Graphs:
                
                if Graph in list( Label_Dictionary.keys() ):
                    
                    plt.figure( dpi = 200 )
                    
                    plt.ylabel( Label_Dictionary[ Graph ] )

                    # TOP CELL FIXED CASE:
                        
                    if Parameters[ "Fixed Multi-Junction Cell" ] == "Top":
                        
                        if Parameters[ "Bottom Cell Material Model" ] == "Simulated":

                            plt.plot( BC_E_gs, Output[ Graph ] )  
                        
                        else:
                            
                            plt.plot( BC_E_gs.values(), Output[ Graph ], '.' )  
                        
                        plt.xlabel( 'Bottom Cell Gap, $E_{\\mathrm{g},2}$ (eV)' )
                        
                    # BOTTOM CELL FIXED CASE:

                    if Parameters[ "Fixed Multi-Junction Cell" ] == "Bottom":
                        
                        if Parameters[ "Top Cell Material Model" ] == "Simulated":

                            plt.plot( E_gs, Output[ Graph ] )  
                        
                        else:
                            
                            plt.plot( E_gs.values(), Output[ Graph ], '.' )  

                        plt.xlabel( 'Top Cell Gap, $E_{\\mathrm{g},1}$ (eV)' )
                    
                    if Graph == 'EQE_EL':
                        
                        plt.yscale( 'log' ) 
                    
                    plt.savefig( path.join( Output_Path, Graph + '.jpg' ) )

    print( "Simulation Complete." )
    
Label_Dictionary = { 'E_g' : '$E_{\\mathrm{g}}$ (eV)',
                    
                     'E_lower' : '$E_{\\mathrm{lower}}$ (eV)',
                    
                    'Intensity [mW/cm2]' : 'Intensity (mW cm$^{-2}$)',
                    
                    'J_sc' : '$J_{\\mathrm{sc}}$ (mA cm$^{-2}$)',
                    
                    'V_oc' : '$V_{\\mathrm{oc}}$ (V)',
                    
                    'V_oc_rad' : '$V_{\\mathrm{oc}}^{\\mathrm{rad}}$ (V)',

                    'V_mpp' : '$V_{\\mathrm{mpp}}$ (V)',

                    'J_mpp' : '$J_{\\mathrm{mpp}}$ (mA cm$^{-2}$)',

                    'P_mpp' : '$P_{\\mathrm{mpp}}$ (mW cm$^{-2}$)',
                    
                    'FF' : 'FF',
                    
                    'PCE' : 'PCE',
                    
                    'EQE_EL' : '$\\mathrm{EQE}_{\\mathrm{EL}}$',
                    
                    'Delta_V_oc_nr' : '$\\Delta V_{\\mathrm{oc}}^{\\mathrm{nr}}$ (V)',

                    'Area_cell' : '$A_{\\mathrm{cell}}$ (cm$^2$)',
                    
                    'I_sc_Module' : '$I_{\\mathrm{sc, Module}}$ (mA)',

                    'I_mpp_Module' : '$I_{\\mathrm{mpp, Module}}$ (mA)',

                    'V_oc_Module' : '$V_{\\mathrm{oc, Module}}$ (V)',

                    'V_mpp_Module' : '$V_{\\mathrm{mpp, Module}}$ (V)',

                    'P_mpp_Module' : '$P_{\\mathrm{mpp, Module}}$ (mW)',

                    'PCE_Module' : '$\\mathrm{PCE}_{\\mathrm{Module}}$',

                    'eta_interp' : '$\\eta_{\\mathrm{EHB, interp}}$',

                    'eta_system' : '$\\eta_{\\mathrm{system}}$'

                    }          
      
#==============================================================================
# If GUI version, call above function when clicked; else, execute it
#==============================================================================

if GUI:    
     
    def Run_Sims_on_Click( Button ):
        
        Value_Storer( '' )
        
        Simulation_Controller()
        
        GUI_Graph_Generator( '' )
        
        Copy_Output.disabled = False
    
    Run_Simulations_Button.on_click( Run_Sims_on_Click )

else:
    
    Simulation_Controller()