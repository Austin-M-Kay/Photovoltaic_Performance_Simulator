#!/usr/bin/env python
# coding: utf-8

# A Tool for Estimating Photovoltaic Performance Under Arbitrary Illumination Conditions

# Austin M. Kay, Maura Fitzsimons, Gregory Burwell, Paul Meredith, Ardalan Armin, & Oskar, J. Sandberg

# Sustainable Advanced Materials (Sêr-SAM), Department of Physics, Swansea University, Singleton Park, Swansea SA2 8PP, United Kingdom

# __Email__: (A. M. Kay) 954708@swansea.ac.uk; o.j.sandberg@swansea.ac.uk; ardalan.armin@swansea.ac.uk;

# This computational tool was created to estimate the performance of a photovoltaic devices under arbitrary illumination conditions. The default spectra that the tool has been made with are the 'warm white' 2700K LED, the 'cool white' 4000K LED, the CIE Standard Illuminant LED-B4, and the standard AM1.5G spectrum for sunlight. Additional specta can be added as new sheets in the "Spectra.xlsx" Excel file. The tool can be used to estimate performance using either a simulated photovoltaic external quantum efficiency (EQE_PV) spectrum (including step functions, sub-gap Urbach tails, and detailed models for organic semiconductor absorption), or an experimentally-determined EQE_PV spectrum (which must be saved in an Excel file in the "EQE_Spectra" folder); the latter will give a realistic estimate for the photovoltaic performance of real devices.

# Two types of simulation can be performed with this tool: (i) figures-of-merit versus the optical gap, and (ii) figures-of-merit versus the intensity of the incident light. For both simulated and experimental EQE_PV, non-radiative open-circuit voltage losses can be included in the simulations using different models, including an emprical model contributed by Maura Fitzsimons. These losses are explored further in the manuscript this computational tool accompanies.

# The following content may be broken into four main sections, followed by the user interface. Firstly, the theoretical background supporting the tool is outlined (alongside the functions that implement it) in Section 1. Following this, in Section 2, the irradiance spectra are loaded into the script and converted to photon flux spectra. Next, in Section 3, the widgets needed to create the User Interface are defined, linked, and compiled. Following this, in Section 4, additional tools used to support the code's operation are defined, including an interpolating function. Finally, in Section 5, the user interface is presented; to quickly load and use the user interface, select "Cell -> Run All" from the top of the screen.

# While every care has been taken to ensure no faults or bugs are present in this tool, we would like to correct any issues users may have - please report these to (A.M. Kay) 954708@swansea.ac.uk. The most up-to-date version of this tool (and its manual) will always be available at https://github.com/Austin-M-Kay).

# The following widget is defined to keep track of script-initialising progress:

from ipywidgets import *

Progress = IntProgress( value = 0 ,
                       
                       min = 0 ,
                       
                       max = 100 , 
                       
                       description = 'Initialisation Progress: ',
                      
                       style = { 'description_width' : 'initial' } )

display( Progress )

# 1. Theoretical Background

# In the following section, the theoretical background is interwoven with accompanying Python functions, ultimately culminating in the calculation of the power conversion efficiency, PCE. 

# 1.1. Short-Circuit Current Density

# To begin, the short-circuit current density, J_sc is defined as [[1](#Ref_Nelson),[2](#Ref_Wurfel)]

# J_sc = q\int_0^\infty EQE_PV(E)\cdot\it{\Phi}_\mathrm{source}(E)\,\mathrm{d}E,

# where q is the elementary charge, EQE_PV(E) is the photovoltaic external quantum efficiency (the ratio of free charge carriers out to the number of photons in) at a given photon energy E, and \it{\Phi}_\mathrm{source}(E) is the spectral photon flux of the incident light. To implement this equation into the code, the following Python function is defined for evaluating the short-circuit current density. This function utilises an in-built funtion from the SciPy library, "simps", which evaluates an integral numerically using Simpson's rule.

from scipy.integrate import simps

# Furthermore, the elmentary charge of the electron q is imported from the SciPy library using:

from scipy.constants import e

# These components are combined to define the following function for evaluating the short-circuit current density - the input parameters for this function are a series of wavelengths, the corresponding EQE_PV values, and the corresponding photon fluxes. The latter two are expected to be provided at the same wavelengths, whether that be measured or interpolated. 

def Short_Circuit_Current_Density_Calaculator( Wavelengths , EQEs , Photon_Fluxes ):

    """This function calculates the short-circuit current density using an EQE_PV spectrum (unitless), its wavelengths (in
    
    nm), and the spectral photon flux of the incident light (in units of mW/cm2/nm). The dimension of each input should be 
    
    the same; i.e., the photon flux and EQE_PV should be provided at the same wavelengths. These quantities are also
    
    expected to be stored in NumPy arrays, such that their product may readily be taken. If this is not the case, and lists 
    
    or tuples are give, this function will convert the lists or tuples to NumPy arrays."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Data Type of Input Data
    #-----------------------------------------------------------------------------------------------------------------------
    
    Wavelengths_Type = type( Wavelengths ).__name__ 
    
    EQE_Type = type( EQEs ).__name__ 
    
    Flux_Type = type( Photon_Fluxes ).__name__ 
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Correct Data Type If Need Be
    #-----------------------------------------------------------------------------------------------------------------------

    if Wavelengths_Type == 'list' or EQE_Type == 'list' or Flux_Type == 'list':
           
        Wavelengths = array( Wavelengths )
        
        EQEs = array( EQEs )
        
        Photon_Fluxes = array( Photon_Fluxes )
            
    if Wavelengths_Type == 'tuple' or EQE_Type == 'tuple' or Flux_Type == 'tuple':
           
        Wavelengths = array( Wavelengths )
        
        EQEs = array( EQEs )
        
        Photon_Fluxes = array( Photon_Fluxes )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Evaluate Short-Circuit Current Density and Return Value
    #-----------------------------------------------------------------------------------------------------------------------
    
    return e * simps( y = EQEs * Photon_Fluxes , x = Wavelengths )

# 1.2. Dark Saturation Current Density

# The dark saturation current density (in the radiative limit), J_0^rad, holds a very similar definition to the short-circuit current density,

# J_0^rad = q\int_0^\infty EQE_PV(E)\cdot\it{\Phi}_bb(E)\,\mathrm{d}E,


# where the Planck black-body radiation spectrum (i.e., the number of unpolarised photons radiated per unit area per unit time per energy interval into a hemispherical solid angle) is given by

# \it{\Phi_bb}(E) = \frac{2\pi E^2}{h^3c^2}[\exp(\frac{E}{k_BT}\right) - 1\right]^{-1}.


# where h is the Planck constant, k_B is the Boltzmann constant, T is the temperature, and the Bose-Einstein distribution has been approximated as the Maxwell-Boltzmann approximation to obtain the right-hand side (valid for E>3k_BT). As the spectra are usually provided in terms of photon wavelength rather than energy, it is useful to instead define the spectral photon flux density as

# \it{\Phi}_\mathrm{Planck}(\lambda)\,\mathrm{d}\lambda = \Phi_\mathrm{Planck}(E)\,\mathrm{d}E=\Phi_\mathrm{Planck}(E)\,|\frac{\mathrm{d}E}{\mathrm{d}\lambda}\right|\mathrm{d}\lambda.

# Then, as E=\frac{hc}{\lambda}\Rightarrow |d_\lambda E|=\frac{hc}{\lambda^2}, Equation ([4](#Planck2)) gives

# \it{\Phi}_\mathrm{Planck}(\lambda) =\frac{2\pi c}{\lambda^4}\frac{1}{\exp(\frac{hc}{\lambda k_BT}\right)-1}.

# Both forms of the Planck spectrum are implemented into the code through the following functions:

from numpy import *
from scipy.constants import e, h, c, k

def Planck_Photon_Flux_Energy( Photon_Energy , Temperature ):
    
    """At a given photon energy (in eV) and temperature (in K), determine the photon flux as defined by the Planck spectrum
    
    in the limit that the Maxwell-Boltzmann approximation is valid. Give this flux in units of #photons/cm2/s/ev."""
    
    Pre_Factor = 2 * pi * ( e * Photon_Energy ) ** 2 / h ** 3 / c ** 2  # Units are #photons/m2/s/J
    
    Pre_Factor = Pre_Factor * 1e-4 / e # Convert to units of #photons/cm2/s/eV
    
    return 1e3 * Pre_Factor / ( exp( e * Photon_Energy / k / Temperature ) - 1 )

def Planck_Photon_Flux_Wavelength( Wavelength , Temperature ):
    
    """At a given photon wavelength (in nm) and temperature (in K), determine the photon flux as defined by the Planck 
    
    spectrum in the limit that the Maxwell-Boltzmann approximation is valid. Give this flux in units of 
    
    #photons/cm2/s/nm."""
    
    Pre_Factor = 2 * pi * c / (1e-9 * Wavelength ) ** 4 * 1e-4 * 1e-9
        
    return 1e3 * Pre_Factor / ( exp( h * c / ( 1e-9 * Wavelength * k * Temperature  ) ) - 1 )

# Using these functions for the Planck black-body flux, the photon flux can be determined at a given temperature and wavelength/energy. Using these values, the dark saturation current in the radiative limit can be determined using the following function:

def Dark_Saturation_Circuit_Current_Density_Calaculator_Rad( Wavelengths , EQE_Spectrum , Photon_Flux_Spectrum ):

    """Compute the dark saturation current density (in mA/cm2) using the wavelengths (in nm) that the input EQE_PV spectrum 
    
    (unitless) and with the input photon flux spectrum (in units of mW/cm2/nm) are given at. The EQE spectrum and the photon
    
    flux spectrum are expected to already be interpolated, such that the values are given at the same points. They are also
    
    expected to be stored in NumPy arrays, such that their product may readily be taken."""
        
    # SciPy's in-built function "simps" has been used to numerically integrating x and y data using Simpson's Rule
    
    return e * simps( y = EQE_Spectrum * Photon_Flux_Spectrum , x = Wavelengths )

# The dark saturation in the non-radiative limit is given by J_0 = \frac{J_0^rad}{EQE_EL}, where the electroluminscent external quantum efficiency, EQE_EL, relates to the non-radiative open-circuit voltage loss, \Delta V_oc^nr, via EQE_EL=\exp(-\frac{q\Delta V_oc^nr}{k_BT}\right). The dark saturation current is evaluated using the following function:

def Dark_Saturation_Circuit_Current_Density_Calaculator( J_0_rad , Delta_V_oc_nr , Temperature ):

    """Determine the dark saturation current density using dark saturation current in the radiaitive limit, the temperature,
    
    and the non-radiative open circuit voltage loss."""
            
    return J_0_rad * exp( e * Delta_V_oc_nr / k / Temperature )

# 1.3. Open-Circuit Voltage

# To express the open circuit voltage V_oc in terms of the current densities, the ideal diode equation is first assumed (with ideal factor m=1). This gives the total current density for a solar cell under illumination J_\mathrm{light} as [[1](#Ref_Nelson)]

# J_\mathrm{light}(V_app)=J_0[\exp(\frac{qV_app}{k_BT}\right)-1\right]-J_sc.

# Here, V_app is the applied voltage, the short circuit current density is defined in Equation ([1](#J_sc_equation)), and the dark saturation current is defined in Equation ([2](#J_0_equation)). A function for determining the total current density using this equation is defined as:

def Total_Current_Density( Applied_Voltage , Dark_Saturation_Current , Short_Circuit_Current , Temperature ):
    
    """Calculate the total current density using the dark saturation current density, the short-circuit current density,
    
    the open-circuit voltage, the temperature, and the constants e and k, being the elementary charge and the Boltzmann 
    
    constant, resepectively."""
    
    return Dark_Saturation_Current * ( exp( e * Applied_Voltage / k / Temperature ) - 1 ) - Short_Circuit_Current

# Under open-circuit conditions where V_app=V_oc, there is no net current generated by the diode (J_\mathrm{light}=0). Consequently,

# V_oc=\frac{k_BT}{q}\ln[\frac{J_sc}{J_0}+1\right].

# The following function is used by the script to calculate the open-circuit voltage:

def Open_Circuit_Voltage( Short_Circuit_Current , Dark_Saturation_Current , Temperature ):
    
    """Determine the open-circuit voltage under the assumption of an ideal diode. The units of the input current densities
    
    aren't important, they just need to be the same such that their ratio is unitless. The temperature should be in 
    
    Kelvin."""
    
    return k * Temperature / e * log( Short_Circuit_Current / Dark_Saturation_Current + 1 )

# As stated above, non-radiative open-circuit voltage losses (\Delta V_oc^nr) may be accounted for using the electroluminescent quantum efficiency

# V_oc=\frac{k_BT}{q}\ln[\frac{J_sc}{J_0^rad}EQE_{EL}+1\right], 

# where EQE_EL=\exp(-\frac{q\Delta V_oc^nr}{k_BT}\right).

# 1.4. Fill Factor

# The fill factor FF is defined as the ratio of the power generated at the maximum power point to the maximum power that could be generated [[2](#Ref_Wurfel)],

# FF=\frac{J_mppV_mpp}{J_scV_oc}.
 
# Here J_mpp is the current density at the maximum power point (mpp), where the applied voltage is V_mpp. This maximum power point is situated at the stationary point that satisfies

# .\frac{\partial P}{\partial V_app}\right|_{V=V_mpp} = .\frac{\partial(J_\mathrm{light} V_app) }{\partial V_app}\right|_{V=V_mpp} = J_\mathrm{light}(V_mpp) + \frac{qV_mpp}{k_BT} J_0 \exp( \frac{qV_mpp}{k_BT}\right) = 0.

# Collecting terms and rearranging gives

# (1+\frac{qV_mpp}{k_BT}\right)\exp( \frac{qV_mpp}{k_BT}\right)=\frac{J_sc}{J_0}+1.

# Making use of the definition of the open-circuit voltage gives:

# (1+\frac{qV_mpp}{k_BT}\right)\exp(1+\frac{qV_mpp}{k_BT}\right)=\exp(1+\frac{qV_oc}{k_BT}\right).


# This equation can be solved (numerically) using the Lambert W function [[3](#Ref_Valluri)], x=W[y], which is defined as the solution to the equation y=xe^x. Using this gives

# V_mpp=\frac{k_BT}{q}\{W[\exp(1+\frac{qV_oc}{k_BT}\right)\right]-1\right\}.

# This equation is implemented into the code using the following function:

from scipy.special import lambertw as W     # Preliminary import of the Lambert W function

def Maximum_Power_Point_Voltage( Open_Circuit_Voltage , Temperature ):
    
    """Determine the voltage that gives the maximum power using the open circuit voltage V_oc in volts and the temperature
    
    in Kelvin. This voltage is evaluated using the Lambwert W function."""
    
    kT = k * Temperature / e
    
    return kT * ( W( exp( 1 + Open_Circuit_Voltage / kT ) ) - 1 ).real # Want real number (shouldn't be complex part anyway)

# With the voltage that gives the maximum power point defined, the current density at this point can be calculated using Equation ([6](#ideal)) with V_app=V_mpp. A function for calculating J_mpp is defined below:

def Maximum_Power_Point_Current_Density( Maximum_Power_Point_Voltage, Dark_Saturation_Current, Short_Circuit_Current, 
                                        
                                        Temperature ):
    
    """Determine the current density at the maximum power point in the units of J_0 and J_sc (should be the same) at the 
    
    voltage that gives the maximum power point V_mpp (in V)."""
    
    return Total_Current_Density( Maximum_Power_Point_Voltage, Dark_Saturation_Current, Short_Circuit_Current, 
                                 
                                 Temperature )

# With the voltage and current at the maximum power point now determined, the maximum power output can be computed using P_mpp=J_mppV_mpp. This is encoded in its own function:

def Maximum_Power_Output( Maximum_Power_Point_Voltage , Maximum_Power_Point_Current_Density ):
    
    """Determine the absolute value of power output at the maximum power point using the voltage and current at this 
    
    point. This quantity will have the units of J_mpp * V_mpp, e.g., mW/cm2."""
    
    return abs( Maximum_Power_Point_Voltage * Maximum_Power_Point_Current_Density )

# The fill factor FF is then given by the ratio of the power generated at the maximum power point to the product of the short-circuit current density and the open-circuit voltage,

# FF = \frac{P_mpp}{J_scV_oc}.

# Equation ([13](#FF_equa)) is implemented into the code using the following function:

def Fill_Factor( Maximum_Power_Point_Current_Density , Short_Circuit_Current_Density , Open_Circuit_Voltage ):
    
    """Calculate the fill factor using the maximum power point power, the short-circuit current density, and the open 
    
    circuit voltage. The units of the power at the maximum power point are expected to equal those of the product of  the
    
    short-circuit current density and the open circuit votlage, such that the fill factor is a unitless quantity with a 
    
    value between nought and one."""
    
    return abs( Maximum_Power_Point_Current_Density / Short_Circuit_Current_Density / Open_Circuit_Voltage )

# 1.5. Power Conversion Efficiency

# The power conversion efficiency \mathrm{PCE} is given by

# \mathrm{PCE} = \frac{P_mpp}{P_\mathrm{source}},

# where the power of the incident light per unit area, P_\mathrm{source}, is given by the integrated spectral irradiance,

# P_\mathrm{light}=\int_0^\infty E\it{\Phi}_\mathrm{source}(E)\,\mathrm{d}E.

# Equations ([14](#PCE_equation)) and ([15](#P_light)) are implemented into the code using the following functions:

def Power_Conversion_Efficiency( Maximum_Power_Point_Current_Density , Light_Power ):
    
    """Calculate the power conversion efficiency using the maximum power point power and the power of the incident light. 
    
    The units of the power at the maximum power point are expected to equal those of the power of the light, such that the 
    
    power conversion efficicency is a unitless quantity with a value between nought and one."""
    
    return abs( Maximum_Power_Point_Current_Density / Light_Power )

def Light_Power( Wavelengths , Irradiances ):
    
    """Determine the total light power of an irradiance spectrum."""
    
    return simps( Irradiances , x = Wavelengths )

# 1.6. Lux Calculation

# For an irradiance spectrum I_\mathrm{source}(E)=E\it{\Phi}_{\mathrm{source}}, with units of power per unit energy per energy interval (or wavelength interval), the illuminance L_\mathrm{source} (in units of lumens per unit area) is given by [[4](#Ref_Sze)]

# L_\mathrm{source}=L_0\int_0^\infty V(E)I_\mathrm{source}(E)\,\mathrm{d}E,

# where V(E) is the luminous efficiency at a given photon energy and L_0 = 683\,\mathrm{lm}\cdot\mathrm{W}^{-1} is a constant; here \mathrm{lm} is the lumen, the unit of luminous flux. By normalising an irradiance spectrum to its total power P_\mathrm{light} as defined in Equation ([15](#P_light)), then a spectrum's value in lux can be related to its value in power per unit aread via

# L_\mathrm{source}=\mathcal{K}P_\mathrm{source},

# where the constant of proportionality \mathcal{K}= L_0\int_0^\infty V(E)I_\mathrm{source}(E)\,\mathrm{d}E holds a unique value for each irradiance spectrum. Once determined, scaling lux values is becomes simple as the constant of proportionality between a given lux value and the corresponing irradiance value is known. Equations for determining lux constants and lux values are defined below: 

def Lux_Constant( Wavelengths , Irradiances , Luminous_Efficiencies ):
    
    """For a given irradiance spectrum, compute the constant of proportionality using the luminous efficiencies (which 
    
    are assumed to be at the same wavelengths). All quantities should be stored in arrays, such that their products may be
    
    readily computed."""
    
    return 683 * simps( y = Luminous_Efficiencies * Irradiances , x = Wavelengths )
    
def Lux_Value( Constant_of_Proportionality , Total_Irradiance_Value ):
    
    """Determine a spectrum's lux using the constant of proportionality for a given spectrum, and its total irradiance."""
    
    return Constant_of_Proportionality * Total_Irradiance_Value


# Update initialisation progress widget:

Progress.value = 5 

# 1.7. Bibliography

# [1] Nelson, J.A., The Physics of Solar Cells. 2003: World Scientific Publishing Company.

# [2] Würfel, P. and U. Würfel, Physics of Solar Cells: From Basic Principles to Advanced Concepts. 2016: John Wiley & Sons.

# [3] Valluri, S.R., D.J. Jeffrey, and R.M. Corless, Some Applications of the Lambert W Function to Physics. Canadian Journal of Physics, 2000. 78(9): p. 823-831.

# [4] Sze, S.M., Y. Li, and K.K. Ng, Physics of Semiconductor Devices. Fourth Edition ed. 2021: John Wiley & Sons, Inc.

# 2. Loading in Spectra

# In this section, the necessary spectra are loaded in. To begin, the tools for loading data are imported:

from os import getcwd, path, listdir

# Following this, the current working directory in which this script is saved is identified the "getcwd" function:

Current_Working_Directory = getcwd()

# The spectra should be stored in the "Spectra" folder of the current directory. The path to this folder is created using by appending the string "Spectra" to the current working directory:

Spectra_Folder_Path = path.join( Current_Working_Directory , "Spectra" )

# The content of the spectra folder is then identified using the function "listdir". This will produce a list of all the files saved in that folder, where each file is callable using its string:

Spectra_Folder_Contents = listdir( Spectra_Folder_Path )

# For each of the files in the list, the path needed to reach that file (which will be useful later) is now defined. These paths are stored in a Python dictionary wherein each element is stored according to a string (in this case, the file's name): 

File_Path_Dictionary = { File_Name :
                       
    path.join( Spectra_Folder_Path , File_Name )
                       
    for File_Name in Spectra_Folder_Contents }

# The file containing the AM1.5G, warm LED, and cool LED spectra is defined to have the following name (which may be modified by user):

AM_LED_Spectra_File_Name = 'Spectra.xlsx'

# With the above now defined, it is removed from the list of all files to leave any additional spectra the user may have added.

Spectra_Folder_Contents.remove( AM_LED_Spectra_File_Name )

# 2.1. Creating a Pandas Data Frame

# A Python module, Pandas, is now used to create a "data frame" - which can be used to select data from a large array using row number and/or column title. To begin with, the necessary modules are loaded:

from pandas import ExcelFile, read_excel

# Pandas' "ExcelFile" function is now applied to load-in the spectra and store them in a data frame. To do this, the file path defined in the previous section is now called from the file path dictionary:   

AM_LED_Data = ExcelFile( File_Path_Dictionary[ AM_LED_Spectra_File_Name ] )

# Following this, the Excel file's sheet names (i.e., all the spectra types saved in the file) are determined using:

AM_LED_Sheet_Names = AM_LED_Data.sheet_names 

# For each of the sheet names in the Excel file, the data is stored in an individual data frame - these data frames are then stored in a dictionary (where each spectrum's data is stored according to its name). The first row in each sheet is to be skipped (hence 'skiprows = 1') and each column has two headers (name and unit - 'header = [0, 1]').

Data_Frames_Dictionary = { Sheet_Name : 
                         
    AM_LED_Data.parse( 
       
        Sheet_Name,
                     
        skiprows = 1, 
                     
        header = [ 0 , 1 ] ) 
                         
    for Sheet_Name in AM_LED_Sheet_Names }

# With the data now loaded in, it can be used. Firstly, the minimum and maximum wavelengths are determined and stored in a dictionary:

Minimum_Spectra_Wavelengths = { Sheet_Name :
                              
    min( Data_Frames_Dictionary[ Sheet_Name ].loc[ : ,  'Wavelength' ].values )
                              
    for Sheet_Name in AM_LED_Sheet_Names }

Maximum_Spectra_Wavelengths = { Sheet_Name :
                              
    max( Data_Frames_Dictionary[ Sheet_Name ].loc[ : ,  'Wavelength' ].values )
                              
    for Sheet_Name in AM_LED_Sheet_Names }

# Update initialisation progress widget

Progress.value = 10 

# In the final part of this section, the irradiance spectra are plotted. Note that these spectra are plotted as they come; i.e., there is currently no re-scaling to plot them clearly together:

import matplotlib.pyplot as plt

# 2.2. Determining Photon Flux

# To determine the photon flux corresponding to a particular irradiance value, the following function is defined to determine the photon energy at a given wavelength:

def Energy_Wavelength_Converter( Energy_or_Wavelength ):
    
    """Convert from photon energy to wavelength or vice versa"""
    
    return h * c / ( e * Energy_or_Wavelength * 1e-9 )

# Following this, for each sheet of the Excel file, the photon fluxes and irradiances are determined, then stored in a dictionary which will be used by the script to determine figures-of-merit.

Photon_Irradiance_Spectra = {}

Photon_Flux_Spectra = {}

for Sheet_Name in AM_LED_Sheet_Names:
    
    Wavelengths = [ Value[ 0 ] for Value in Data_Frames_Dictionary[ Sheet_Name ].loc[ : ,  'Wavelength' ].values ]
    
    Irradiances = [ Value[ 0 ] for Value in Data_Frames_Dictionary[ Sheet_Name ].loc[ : ,  'Intensity per wavelength' ].values ]
    
    Fluxes = [ 0.1 * Irradiances[ i ] / e / Energy_Wavelength_Converter( Wavelengths[ i ] ) for i in range( len( Wavelengths ) ) ]
          
    Photon_Irradiance_Spectra[ Sheet_Name ] = array( [ Wavelengths , Irradiances ] )
    
    Photon_Flux_Spectra[ Sheet_Name ] = array( [ Wavelengths , Fluxes ] )

# 2.3. Determining Total Light Power and Lux

# The powers densities associated with these spectra are determined using Equation ([15](#P_light)), where all the imported irradiance spectra have units of \mathrm{W}\cdot\mathrm{m}^{-2}\cdot\mathrm{nm}^{-1}. Integrating these irradiance spectra therefore gives a total incident power in units of \mathrm{W}\cdot\mathrm{m}^{-2} and so, an addional scale factor of 10^{-1} is included to convert to units of mW/cm2 in the following calculations:

P_lights = { Sheet_Name : 
            
    Light_Power( Photon_Irradiance_Spectra[ Sheet_Name ][ 0 ] , Photon_Irradiance_Spectra[ Sheet_Name ][ 1 ] ) / 10 
    
    for Sheet_Name in AM_LED_Sheet_Names }

# Each of the photon flux spectra is then normalised to its power density and stored in a dictionary using:

Normalised_Photon_Flux_Spectra = {}

for Sheet_Name in AM_LED_Sheet_Names:
    
    Wavelengths , Fluxes = Photon_Flux_Spectra[ Sheet_Name ]
    
    Normalised_Photon_Flux_Spectra[ Sheet_Name ] = array( [ Wavelengths , Fluxes / P_lights[ Sheet_Name ] ] )

# The luminous efficiency data is loaded in using the following code:

Luminous_Efficiency_File_Name = 'Luminous_Efficiency_Data.xlsx'

Luminous_Efficiency_Data = ExcelFile( 
    
    path.join( 
        
        Current_Working_Directory, 
        
        Luminous_Efficiency_File_Name ) )

# This file contains luminous efficiency spectra for each of the three light cones that form the human eye: L - long wavelengths (red), M - medium wavelengths (green), and S - short wavelengths (blue). The available spectra (inferred from the sheet names) are given user-friendly names using the following dictionary:

V_Names = { '2_Deg_V' : 'V 2-deg',
           
           '10_Deg_V' : 'V 10-deg',
           
           'Light_Cone_L' : 'L',
          
          'Light_Cone_M' : 'M',
          
          'Light_Cone_S' : 'S'}

# The spectra are then parsed into individual dictionaries using: 

Luminous_Efficiency_Dictionary = { V_Names[ Sheet_Name ] : 
                         
    Luminous_Efficiency_Data.parse( 
       
        Sheet_Name ) 
                         
    for Sheet_Name in Luminous_Efficiency_Data.sheet_names }

# For each of these data frames, invalid numbers are removed using:

for Key in list( Luminous_Efficiency_Dictionary.keys() ):
    
    Luminous_Efficiency_Dictionary[ Key ] = Luminous_Efficiency_Dictionary[ Key ].dropna()

# Progress update:

Progress.value = 20 

# Each of these dictionaries contain a finite amount of data - extrapolation will be needed to estimate the luminous efficiency outside the data set. This extrapolation is modelled as a Gaussian with respect to the photon energy E of the form:

# V(E)= e^{a_0E^2+b_0E+c_0}\propto e^{[A-BE]^2}.

# All parameters are boundless bar the coefficient of the quadratic term, a_0, which is forced to be negative such that the natural logarithm of the luminous efficiency is an inverse parabola. The parameter c_0 will relate to the amplitude the Gaussian will take at E=0. In terms of the photon wavelength, this model takes the form

# V(E)= \exp(\frac{a_0h^2c^2}{\lambda^2}+\frac{b_0hc}{\lambda}+c_0\right),

# where a_0h^2c^2 must have units of length-squared, b_0hc must have units of length, and c_0 must be unitless.

# To perform the curve fitting, the curve fit function is imported from the SciPy module:

from scipy.optimize import curve_fit 

# To fit the data, the natural logarithm will be taken, such that the data will be fit with a quadratic model of the form:

def Quadratic( x , a_0 , b_0 , c_0 ):
    
    """A model for a quadratic line with coefficients a_0, b_0, and c_0."""
    
    return a_0 * x ** 2 + b_0 * x + c_0 

# The User can change the number of data points that should be fit with the model below:

Number_of_Data_Points_to_Fit = 30

# The parameters used in these Gaussian fittings are to be stored in the following dictionary:

Gaussian_Fit_Parameter_Dictionary = {}

# Using the first and last data points (default thirty - user can modify), the natural logarithm of the luminous efficiency spectra are fit with the above model using the following function:

def Luminous_Efficiency_Extrapolater( Number_of_Data_Points_to_Fit , Wavelengths , Luminous_Efficiencies ):
    
    """For a given data set, extraploate the first and last portion of the luminous efficiencies using the wavelengths. The
    
    wavelengths should be in units of nm and ascending."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Prerequisite Data Sorting
    #-----------------------------------------------------------------------------------------------------------------------
    
    Energies = h * c / ( e * 1e-9 * Wavelengths )
    
    Energies_Flipped = Energies[ ::-1 ] # Flip the array so thei
    
    ln_Luminous_Efficiencies = log( Luminous_Efficiencies )
    
    ln_Luminous_Efficiencies_Flipped = ln_Luminous_Efficiencies[ ::-1 ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Fit First Tail
    #-----------------------------------------------------------------------------------------------------------------------
    
    Optimal_Parameters_1 , Optimal_Parameter_Errors_1 = curve_fit( 
    
        Quadratic, 
    
        Energies_Flipped[ :Number_of_Data_Points_to_Fit ],
    
        ln_Luminous_Efficiencies_Flipped[ :Number_of_Data_Points_to_Fit ],
        
        bounds = ( ( -inf ,  -inf , -inf ) , ( 0 , inf , inf ) ) )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Fit Second Tail
    #-----------------------------------------------------------------------------------------------------------------------
    
    Optimal_Parameters_2 , Optimal_Parameter_Errors_2 = curve_fit( 
    
        Quadratic, 
    
        Energies_Flipped[ -Number_of_Data_Points_to_Fit: ],
    
        ln_Luminous_Efficiencies_Flipped[ -Number_of_Data_Points_to_Fit: ],
    
        bounds = ( ( -inf ,  -inf , -inf ) , ( 0 , inf , inf ) ) )
    
    Scale_Factors = [ ( h * c / ( e * 1e-9 ) ) ** 2 , ( h * c / ( e * 1e-9 ) ) , 1 ]
    
    return Scale_Factors * Optimal_Parameters_1 , Scale_Factors * Optimal_Parameters_2

# For each of the available spectra, the tails are now extrapolated and plotted using: <br/><br/>

Colours = [ 'k' , 'r' , 'g' , 'b' , 'c' , 'm' , 'y' ]

for Key in list( Luminous_Efficiency_Dictionary.keys() ):
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Specify Curve Colour
    #-----------------------------------------------------------------------------------------------------------------------
    
    Index = list( Luminous_Efficiency_Dictionary.keys() ).index( Key )

    Curve_Colour = Colours[ Index % 7 ]
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Load Data
    #-----------------------------------------------------------------------------------------------------------------------
    
    Data = Luminous_Efficiency_Dictionary[ Key ]
    
    Wavelength_Data = Data[ Data.columns[ 0 ] ]
    
    Luminous_Efficiency_Data = Data[ Data.columns[ 1 ] ]
    
    Lower_Plotting_Wavelengths = linspace( 300 , min( Wavelength_Data ) , 101 )

    Upper_Plotting_Wavelengths = linspace( max( Wavelength_Data ) , 1200 , 201 )

    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Extrapolation Parameters and Store Them
    #-----------------------------------------------------------------------------------------------------------------------
        
    Optimal_Parameters_1 , Optimal_Parameters_2 = Luminous_Efficiency_Extrapolater( 
        
        Number_of_Data_Points_to_Fit,
        
        Wavelength_Data, 
        
        Luminous_Efficiency_Data )
    
    Gaussian_Fit_Parameter_Dictionary[ Key ] = [ Optimal_Parameters_1 , Optimal_Parameters_2 ]
    
    # Update Progress
    
    Progress.value = 20 + int( 20 * ( list( Luminous_Efficiency_Dictionary.keys() ).index( Key ) + 1 ) / len(
        
        Luminous_Efficiency_Dictionary.keys() ) )   
    
# In the above plot, the dashed lines indicate the imported data, whereas the solid lines are the extrapolated tails. Using these extrapolations, a function is now defined for determining the luminous efficiency at any wavelength:

def Extrapolated_Luminous_Efficiency( Luminous_Efficiency_Type , Wavelength ):
    
    """For a given luminous efficiency extrapolate (or interpolate) the sepctrum to determine its value at a particular 
    
    wavelength."""
    
    LE_Data = Luminous_Efficiency_Dictionary[ Luminous_Efficiency_Type ]
    
    if Wavelength < min( list( LE_Data[ LE_Data.columns[ 0 ] ] ) ):
                                             
        return exp( Quadratic( 1 / Wavelength , *Gaussian_Fit_Parameter_Dictionary[ Luminous_Efficiency_Type ][ 1 ] ) )
                                             
    if Wavelength > max( list( LE_Data[ LE_Data.columns[ 0 ] ] ) ):
                                             
        return exp( Quadratic( 1 / Wavelength , *Gaussian_Fit_Parameter_Dictionary[ Luminous_Efficiency_Type ][ 0 ] ) )
    
    else:
        
        return interp( Wavelength , LE_Data[ LE_Data.columns[ 0 ] ], LE_Data[ LE_Data.columns[ 1 ] ] )

# Using this, the constants of proportionality (for determining lux values) and the lux values themselves are deduced for each of the initial input spectra:

Constants_of_Proportionality = {}

Lux_Values = {}

for Key in AM_LED_Sheet_Names:
    
    Wavelengths = Photon_Irradiance_Spectra[ Key ][ 0 ]
    
    Constants_of_Proportionality_i = {}
    
    Lux_Values_i = {}
    
    for Luminous_Efficiency_Type in list( Luminous_Efficiency_Dictionary.keys() ):
    
        Constants_of_Proportionality_i[ Luminous_Efficiency_Type ] = Lux_Constant(
            
            Wavelengths, 
                     
            Photon_Irradiance_Spectra[ Key ][ 1 ] / ( P_lights[ Key ] * 10 ), 
                     
            [ Extrapolated_Luminous_Efficiency( 
                         
                Luminous_Efficiency_Type, 
                         
                Wavelength ) for Wavelength in Wavelengths ] )     
             
        Lux_Values_i[ Luminous_Efficiency_Type ] = Lux_Value( 
            
            Constants_of_Proportionality_i[ Luminous_Efficiency_Type ], 
            
            P_lights[ Key ] * 10  ) # Multiplied by 10 to get W/m2
        
    Constants_of_Proportionality[ Key ] = Constants_of_Proportionality_i
    
    Lux_Values[ Key ] = Lux_Values_i
    
    # Update Progress
    
    Progress.value = 40 + int( 40 * ( AM_LED_Sheet_Names.index( Key ) + 1 ) / len( AM_LED_Sheet_Names ) )    

# 3. Prequisites for a User Interface

# In this section, some prerequisite widgets needed to create the interface are defined, linked, and compiled. In [Section 3.1](#EQE_Loader), the widgets used to generate the EQE spectrum loading tool are defined. In [Section 3.2](#Lux_Customisation), the widgets needed to select a spectrum and customise its irradiance value are defined, then in [Section 3.3](#Additional_Widgets), any additional widgets are specified.

#  3.1. EQE Spectrum Loading Tool

# As a prerequisite to building the User Interface, all widgets are imported from Jupyter's 'ipywidgets' library:

from ipywidgets import *

# 3.1.1. Determinining EQE Spectra File Paths

# The path to EQE spectra directory is created using the current working directory (defined in the previous section):

EQE_Spectra_Folder_Path = path.join( Current_Working_Directory , "EQE_Spectra" )

# The files in this folder are then identified and stored using:

EQE_Spectra_Folder_Contents = listdir( EQE_Spectra_Folder_Path )

# The corresponding file paths are determined and stored in the following dictionary:

EQE_File_Path_Dictionary = { File_Name :
                       
    path.join( EQE_Spectra_Folder_Path , File_Name )
                       
    for File_Name in EQE_Spectra_Folder_Contents }

# 3.1.2. EQE Spectrum Selection Widgets

# With the paths to the EQE spectra specified, a button widget is now defined for to commence the EQE-spectrum loading process:

Load_EQE_Button = Button( 

    description = 'Add EQE Spectrum' )

# Ultimately, when pressed this button will open a data-loading tool-box such that the experimental spectra can be loaded. Before this can happen, each tool in the toolbox must be created; the first of these is a dropdown widget to select a spectrum from the available files: 

EQE_Choice_Widget_Options = EQE_Spectra_Folder_Contents.copy()

EQE_Choice_Widget = Combobox(

    options = EQE_Choice_Widget_Options,
    
    placeholder = 'File Name',
    
    description = 'Select File Name:',
    
    style = { 'description_width' : 'initial' } )

# In some files, a number of rows may need to be skipped, or the header of the columns may extend across more than one row. These cases are accounted for by defining the following widgets:

Skip_Rows_Widget = BoundedIntText( 

    value = 0,

    min = 0,

    max = 1e40,

    description = 'Number of Rows to Skip:',

    style = { 'description_width' : 'initial' } ) 

Header_Length_Widget = BoundedIntText( 

    value = 1,

    min = 0,

    max = 1e40,

    description = 'Header Length =',

    style = { 'description_width' : 'initial' } ) 

Skip_Rows_Widget.layout.display = 'none' 
        
Header_Length_Widget.layout.display = 'none'

# Following this, the widgets that form the final part of the EQE-selecting box are specified. These include a widget to cancel the data-loading process and a widget to save the EQE spectrum that the User will have selected.

Close_Box_Button = Button( description = 'Cancel' , button_style = 'danger' )

Add_Spectrum_Button = Button( description = 'Add Spectrum' , button_style = 'success' )

Coda_Box = HBox( [ Close_Box_Button , Add_Spectrum_Button ] )

Coda_Box.layout.display = 'none'

# The coda box is then compiled with the EQE-selecting widget into an "Add Spectrum" box, which is set to initially be hidden. When the User clicks the "Add EQE Spectrum" button, the box will be revealed. Note, a blank placeholder label is added into this box such that further tools can take its place in the next section. 

Placeholder_Label = Label()

Error_Label = Valid( 
    
    value = False, 
    
    readout = '',

    layout = Layout(
    
        width = '100%' ) )

Error_Label.layout.display = 'none'

Add_Spectrum_Box = VBox( [ EQE_Choice_Widget , Skip_Rows_Widget , Header_Length_Widget , Placeholder_Label ] )

Add_Spectrum_Box.layout.display = 'none'

Compiled_Add_Spectrum_Box = VBox( [ 
    
    Load_EQE_Button,
    
    Add_Spectrum_Box,
    
    Coda_Box,
    
    Error_Label ])

# A label for displaying error messages was also compiled into the box, which will be useful for error messages later on. The buttons for opening and closing the Load EQE spectrum dialogue box are now tied to operations using the following two functions:

def Open_EQE_Spectrum_Adder( Button ):
    
    """If the user chooses to add a spectrum, open the box."""
    
    Add_Spectrum_Box.layout.display = None
    
    Coda_Box.layout.display = None
    
Load_EQE_Button.on_click( Open_EQE_Spectrum_Adder )

def Close_EQE_Spectrum_Adder( Button ):
    
    """If the user chooses to add a spectrum, open the box."""
    
    Add_Spectrum_Box.layout.display = 'none'
    
    Coda_Box.layout.display = 'none'
    
    Error_Label.layout.display = 'none'
    
Close_Box_Button.on_click( Close_EQE_Spectrum_Adder )

# The save data button's function will be defined at the end of this section.

# 3.1.3. Data Customisation Widgets

# The functions in this section essentially perform one job but each gets a bit involved. We will summarise each of their duties before defining them.

# The following function loads the data from the Excel file the User selects, then proceeds to use other functions to (i) generate a control panel populated with widgets for customising column selection and specifying units, (ii) giving the widgets their instructions such that, e.g., pressing "Update Labels" updates the columns associated with the stored data, and (iii) presents a preview of the loaded data in a Pandas Data Frame:

def On_Change_Load_EQE_Widgets_Updater( Change ):
    
    """If the new filename is valid, generate the widgets and provide a preview of the data frame."""
   
    #-----------------------------------------------------------------------------------------------------------------------
    # Remove Placeholder from 'Add Spectrum' Box 
    #-----------------------------------------------------------------------------------------------------------------------
       
    Add_Spectrum_Box.children = Add_Spectrum_Box.children[ : -1 ] 
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Only Proceed if the User's Choice of EQE Filename Is Valid (i.e., The File is Present in the Folder)
    #-----------------------------------------------------------------------------------------------------------------------
        
    if EQE_Choice_Widget.value in EQE_Choice_Widget_Options:
        
        Skip_Rows_Widget.layout.display = None 
        
        Header_Length_Widget.layout.display = None
        
        #-------------------------------------------------------------------------------------------------------------------
        # Hide the Error Label if it Had Previously Been Revealed
        #-------------------------------------------------------------------------------------------------------------------
        
        Error_Label.readout = ''
        
        Error_Label.layout.display = 'none'
        
        #-------------------------------------------------------------------------------------------------------------------
        # Load the Data into a Globally-Defined Data Frame
        #-------------------------------------------------------------------------------------------------------------------
        
        global Loaded_DataFrame 
        
        if Header_Length_Widget.value != 0:
            
            Header = list( range( Header_Length_Widget.value ) ) 
            
        else:
            
            Header = None
            
        Loaded_DataFrame = read_excel( 
                    
            EQE_File_Path_Dictionary[ EQE_Choice_Widget.value ] ,
                    
            skiprows = Skip_Rows_Widget.value, 
                
            header = Header ) 
        
        #-------------------------------------------------------------------------------------------------------------------
        # Store the Column Titles in a Globally-Defined Variable
        #-------------------------------------------------------------------------------------------------------------------
            
        global Columns
        
        Columns = Loaded_DataFrame.columns
        
        #-------------------------------------------------------------------------------------------------------------------
        # Use Another Function to Generate a Grid of Data-Selecting Widgets
        #-------------------------------------------------------------------------------------------------------------------
    
        global Grid
        
        Grid = Data_Frame_Control_Box_Maker( Loaded_DataFrame )

        #-------------------------------------------------------------------------------------------------------------------
        # Create Widgets to Specify Whether a Variable Has Been Changed (Prevents Crashing)
        #-------------------------------------------------------------------------------------------------------------------
        
        global Variable_Changed_Widgets
        
        Variable_Changed_Widgets = { Column :
                           
            Checkbox( value = False ) 
                           
            for Column in Loaded_DataFrame.columns }    
      
        #-------------------------------------------------------------------------------------------------------------------
        # Generate a Button for Updating the Column Labels, Function Defined Elsewhere
        #-------------------------------------------------------------------------------------------------------------------
        
        global Update_Labels_Button
        
        Update_Labels_Button = Button( description = 'Update Labels' )

        Update_Labels_Button.on_click( Update_Lables_Button_Updater )
        
        #-------------------------------------------------------------------------------------------------------------------
        # Generate the Output DataFrame
        #-------------------------------------------------------------------------------------------------------------------
        
        Output_Table = Output()
        
        Add_Spectrum_Box.children = Add_Spectrum_Box.children + ( Output_Table , )
        
        global Output_DataFrame
        
        Output_DataFrame = Output_DataFrame_Generator( Loaded_DataFrame , Column_Labels , Include_Checkboxes )
        
        with Output_Table:
            
            display( Update_Labels_Button , Grid , Output_DataFrame )    
                
        for Column in Loaded_DataFrame.columns:
    
            Include_Checkboxes[ Column ].observe( Column_Unincluder , names = 'value' )
    
            Independent_Variable_Checkboxes[ Column ].observe( Independent_Variable_Changer , names = 'value' )
    
            Dependent_Variable_Checkboxes[ Column ].observe( Dependent_Variable_Changer , names = 'value' )
                
    else:
        
        Add_Spectrum_Box.children = Add_Spectrum_Box.children + ( Placeholder_Label , )

# The above function is called when the User changes the value of the EQE file-selecting widget, it does this because of the following instruction:

EQE_Choice_Widget.observe( On_Change_Load_EQE_Widgets_Updater , names = 'value' )

# The function for loading EQE spectra relies on two other functions which, in turn, rely on further functions to give widgets instructions etc. The first of these is more simple, it generates a data frame based on which file the User has selected, their choice of column labels, and which columns have been selected:

def Output_DataFrame_Generator( Loaded_DataFrame , Column_Labels , Include_Checkboxes ):
    
    """Generate the output dataframe."""
    
    Columns = Loaded_DataFrame.columns
    
    Output_DataFrame = Loaded_DataFrame.copy()

    Output_DataFrame.drop( [ Column for Column in Columns if Include_Checkboxes[ Column ].value == False ] , axis = 1 )    
    
    Output_DataFrame.rename( { Column : Column_Labels[ Column ].value for Column in Columns } )
    
    return Output_DataFrame

# The second function is a little more involved. It generates all the widgets that the User can use to select columns from an Excel file, etc. It then compiles them into a grid where each column has a label:

def Data_Frame_Control_Box_Maker( Loaded_DataFrame ):
    
    """Using an input data frame, create all the widgets to control the data loading."""
    
    Columns = Loaded_DataFrame.columns
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Globally Define Some Quantities So They Can Be Referred To Elsewhere
    #-----------------------------------------------------------------------------------------------------------------------
    
    global Column_Labels, Include_Checkboxes, Independent_Variable_Checkboxes, Dependent_Variable_Checkboxes
    
    global Type_Widgets, Unit_Widgets
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Widgets to Specify the Column Labels (Initially, Let This Be the Column Titles of the Excel File)
    #-----------------------------------------------------------------------------------------------------------------------
    
    Column_Labels = { Column :
                    
        Text( 
            
            value = Column,
        
            layout = Layout( 
            
                width = '16.6666666666%' , 
            
                display = 'flex',
            
                justify_content = 'center' ) )
                    
        for Column in Columns }
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Widgets to Include / Remove Widgets     
    #-----------------------------------------------------------------------------------------------------------------------
    
    Include_Checkboxes = { Column : 
                         
        Checkbox(
        
            value = True,
            
            indent = False,
        
            layout = Layout( 
            
                width = '16.6666666666%',
                
                display = 'flex',
            
                justify_content = 'center' ),
        
            style = { '_view_name' : str( Column ) } )
                         
        for Column in Columns }
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Widgets to Set as Independent Variable or Dependent Variable Data     
    #-----------------------------------------------------------------------------------------------------------------------
    
    Independent_Variable_Checkboxes = { Column : 
                         
        Checkbox(
        
            value = False,
        
            indent = False,
            
            layout = Layout( 
            
                width = '16.6666666666%',
                
                display = 'flex',
            
                justify_content = 'center' ),
        
            style = { '_view_name' : str( Column ) } )
                                       
        for Column in Columns }
    
    Dependent_Variable_Checkboxes = { Column : 
                         
        Checkbox(
        
            value = True,
            
            indent = False,
        
            layout = Layout( 
            
                width = '16.6666666666%',
        
                display = 'flex',
            
                justify_content = 'center' ),
        
            style = { '_view_name' : str( Column ) } )
                                     
        for Column in Columns }
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Widgets to Specify Whether an Independent Variable is in Terms of Photon Energy or Wavelength     
    #-----------------------------------------------------------------------------------------------------------------------
    
    Type_Widgets = { Column :
                   
        RadioButtons( 
            
            options = [ 'Wavelength' , 'Energy'],
        
            value = 'Wavelength',
            
            layout = Layout( 
            
                width = '16.6666666666%' , 
            
                display = 'flex',
            
                justify_content = 'center' ), 
        
            disabled = True,
        
            style = { '_view_name' : str( Column ) } )
                    
        for Column in Columns }
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Widgets to Specify the Units of a Dependent Variable   
    #-----------------------------------------------------------------------------------------------------------------------
    
    Unit_Widgets = { Column :
                   
        RadioButtons( 
            
            options = [ 'Unitless' , '%'],
        
            value = 'Unitless',
            
            layout = Layout( 
            
                width = '16.6666666666%' , 
            
                display = 'flex',
            
                justify_content = 'center' ), 
        
            disabled = False,
        
            style = { '_view_name' : str( Column ) } )
                    
        for Column in Columns }
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Assume the First Column Initially Holds the Independent Variable Data
    #-----------------------------------------------------------------------------------------------------------------------
    
    First_Column = Loaded_DataFrame.columns[ 0 ]
    
    Independent_Variable_Checkboxes[ First_Column].value = True
    
    Dependent_Variable_Checkboxes[ First_Column ].value = False
    
    Type_Widgets[ First_Column ].disabled = False
    
    Unit_Widgets[ First_Column ].disabled = True
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create The Column Title Widgets with Centred Text
    #-----------------------------------------------------------------------------------------------------------------------
    
    Column_Titles = HBox( [
        
        Label( 
            
            value = 'Column Label' ,
            
            layout = Layout( 
            
                width = '16.6666666666%' , 
            
                display = 'flex',
            
                justify_content = 'center',
            
                border = 'groove' ) ),
        
        Label( 
            
            value = 'Include Column' , 
            
            layout = Layout( 
                
                width = '16.6666666666%' , 
                
                display = 'flex',
                
                justify_content = 'center',
            
                border = 'groove' ) ),
        
        Label( 
            
            value = 'Independent Variable' , 
            
            layout = Layout( 
                
                width = '16.6666666666%' , 
                
                display = 'flex',
                
                justify_content = 'center',
            
                border = 'groove' ) ),
        
        Label( 
            
            value = 'Type' , 
            
            layout = Layout( 
                
                width = '16.6666666666%' , 
                
                display = 'flex',
                
                justify_content = 'center',
            
                border = 'groove' ) ),
        
        Label(
            
            value = 'Dependent Variable' , 
            
            layout = Layout( 
                
                width = '16.6666666666%' , 
                
                display = 'flex',
                
                justify_content = 'center',
            
                border = 'groove' ) ),
        
        Label( 
            
            value = 'Unit' , 
            
            layout = Layout( 
                
                width = '16.6666666666%' , 
                
                display = 'flex',
                
                justify_content = 'center',
            
                border = 'groove' ) ) ] )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Store All the Output in a Grid
    #-----------------------------------------------------------------------------------------------------------------------
    
    Grid = GridspecLayout( len( Columns ) + 1 , 6 )
    
    Grid[ 0 , : ] = Column_Titles
    
    for i in range( len( Columns ) ):
        
        Column = Columns[ i ]
        
        Grid[ i + 1 , : ] = HBox( [
            
            Column_Labels[ Column ], 
            
            Include_Checkboxes[ Column ],
        
            Independent_Variable_Checkboxes[ Column ],
            
            Type_Widgets[ Column ],
            
            Dependent_Variable_Checkboxes[ Column ],
            
            Unit_Widgets[ Column ] ] )
    
    return Grid

# The following function updates the Pandas Data Frame to remove columns that the User chooses to remove. It also disables all widgets associated with that column except the widget needed to include it again. 

def Column_Unincluder( Change ):
    
    """If the user chooses to not incldue a particular column, disable its customisation options."""
    
    Column = Change.owner.style._view_name
            
    if Change[ 'new' ] == True:
                
        Independent_Variable_Checkboxes[ Column ].disabled = False
    
        Dependent_Variable_Checkboxes[ Column ].disabled = False
        
        if Independent_Variable_Checkboxes[ Column ].value == True:
            
            Type_Widgets[ Column ].disabled = False
            
            Unit_Widgets[ Column ].disabled = True
        
        else:
            
            Type_Widgets[ Column ].disabled = True
            
            Unit_Widgets[ Column ].disabled = False
        
    if Change[ 'new' ] == False:
                
        Independent_Variable_Checkboxes[ Column ].disabled = True
    
        Dependent_Variable_Checkboxes[ Column ].disabled = True
        
        Type_Widgets[ Column ].disabled = True
        
        Unit_Widgets[ Column ].disabled = True
        
#        Output_DataFrame = Output_DataFrame.drop( [ Column_Labels[ Column ].value ] , axis = 1 )
    
    global Output_DataFrame, Columns
    
    Columns = Loaded_DataFrame.columns
    
    Output_DataFrame = Loaded_DataFrame.copy()

    Output_DataFrame = Output_DataFrame.drop( [ Col for Col in Columns if Include_Checkboxes[ Col ].value == False ] , axis = 1 )    
    
    Output_DataFrame = Output_DataFrame.rename( { Col : Column_Labels[ Col ].value for Col in Columns } , axis = 1 )
    
    Add_Spectrum_Box.children = Add_Spectrum_Box.children[ : -1 ] 
        
    Output_Table = Output()
        
    Add_Spectrum_Box.children = Add_Spectrum_Box.children + ( Output_Table , )
    
    with Output_Table:
        
        display( Grid , Output_DataFrame )

# The following function is used to switch between a column between being a dependent variable and an independent variable:

def Independent_Variable_Changer( Change ):
    
    """If the user changes the column from containing independent variable data to dependent variable data, update the
    
    widgets."""
    
    Column = Change.owner.style._view_name
    
    if Variable_Changed_Widgets[ Column ].value != True:
    
        Variable_Changed_Widgets[ Column ].value = True # Prevent the other variable from repeating the process
        
        if Change[ 'new' ] == True:
                    
            Dependent_Variable_Checkboxes[ Column ].value = False
        
            Type_Widgets[ Column ].disabled = False
        
            Unit_Widgets[ Column ].disabled = True     
            
        if Change[ 'new' ] == False:
            
            Variable_Changed_Widgets[ Column ].value = True
        
            Dependent_Variable_Checkboxes[ Column ].value = True
        
            Type_Widgets[ Column ].disabled = True
        
            Unit_Widgets[ Column ].disabled = False  
                    
        Variable_Changed_Widgets[ Column ].value = False
        
    else:
            
        pass

# Whereas the following function does the opposite:

def Dependent_Variable_Changer( Change ):
    
    """If the user changes the column from containing independent variable data to dependent variable data, update the
    
    widgets."""
    
    Column = Change.owner.style._view_name
    
    if Variable_Changed_Widgets[ Column ].value != True:
            
        Variable_Changed_Widgets[ Column ].value = True
        
        if Change[ 'new' ] == True:
        
            Independent_Variable_Checkboxes[ Column ].value = False
        
            Type_Widgets[ Column ].disabled = True
        
            Unit_Widgets[ Column ].disabled = False        
            
        if Change[ 'new' ] == False:
        
            Independent_Variable_Checkboxes[ Column ].value = True
        
            Type_Widgets[ Column ].disabled = False
        
            Unit_Widgets[ Column ].disabled = True
                    
        Variable_Changed_Widgets[ Column ].value = False
    
    else:
            
        pass

# The following function is used to update the labels of the data frame if the user presses the button:

def Update_Lables_Button_Updater( Button ):
    
    """Update the labels of the data frame when clicked."""
    
    if EQE_Choice_Widget.value in EQE_Choice_Widget_Options:    
  
        global Output_DataFrame
    
        Output_DataFrame = Loaded_DataFrame.copy()
        
        Output_DataFrame = Output_DataFrame.drop( [ Col for Col in Columns if Include_Checkboxes[ Col ].value == False ] , axis = 1 )    
    
        Output_DataFrame = Output_DataFrame.rename( { Col : Column_Labels[ Col ].value for Col in Columns } , axis = 1 )
    
        Add_Spectrum_Box.children = Add_Spectrum_Box.children[ : -1 ] 
        
        Output_Table = Output()
        
        Add_Spectrum_Box.children = Add_Spectrum_Box.children + ( Output_Table , )

        with Output_Table:
        
            display( Update_Labels_Button , Grid , Output_DataFrame )        

# 3.1.4. Data Storing Functions

# Once the user has selected their desired EQE, spectra, they will need to press the "Add Spectrum" button. This button will run the following function: 

Count_of_Rows = IntText( value = 0,
                       
                       max = 1e40 )

def Save_Added_EQE_Spectrum( Button ):
    
    """If the user chooses to save the spectrum, do so."""
        
    if EQE_Choice_Widget.value in EQE_Choice_Widget_Options:
        
        #-------------------------------------------------------------------------------------------------------------------
        # If The User Has Selected a Valid Filename, Else Give Error Message
        #-------------------------------------------------------------------------------------------------------------------
        
        if len( Output_DataFrame.columns ) < 2:
            
            #---------------------------------------------------------------------------------------------------------------
            # Ensure Two or More Columns of Data Have Been Provided, Else Give Error Message
            #---------------------------------------------------------------------------------------------------------------
        
            Error_Label.layout.display = None
            
            Error_Label.readout = 'Too Few Columns'
            
        else:
            
            #---------------------------------------------------------------------------------------------------------------
            # Make a Tally of How Many Columns are Independent Variables and How Many are Dependent Variables
            #---------------------------------------------------------------------------------------------------------------
            
            Number_of_Independents = 0
            
            Number_of_Dependents = 0
            
            for Column in Columns:
                
                if Include_Checkboxes[ Column ].value == True:
                
                    if Independent_Variable_Checkboxes[ Column ].value == True:
                    
                        Number_of_Independents += 1
            
                    if Dependent_Variable_Checkboxes[ Column ].value == True:
                    
                        Number_of_Dependents += 1
                        
            #---------------------------------------------------------------------------------------------------------------
            # If None of the Columns are Independent (or Dependent) - Error, Can't Calculate
            #---------------------------------------------------------------------------------------------------------------
        
            if Number_of_Independents == 0 or Number_of_Dependents == 0:
                
                Error_Label.layout.display = None
                
                if Number_of_Independents == 0:
                    
                    Error_Label.readout = 'No Independent Variables'
                
                if Number_of_Dependents == 0:
            
                    Error_Label.readout = 'No Dependent Variables'
                
            else:
                
                #-----------------------------------------------------------------------------------------------------------
                # If Everything is in Order, Proceed - Use the Data Compiler Function to Compile the Data
                #-----------------------------------------------------------------------------------------------------------
                
                global EQE_Spectra_to_Investigate
                
                Compiled_Data = Data_Compiler( Output_DataFrame )
                
                Current_Keys = list( EQE_Spectra_to_Investigate.keys() )
                
                New_Keys = list( Compiled_Data.keys() )
                    
                for i in range( len( New_Keys ) ):
                    
                    Key = New_Keys[ i ]
                    
                    if Key in Current_Keys:
                        
                        New_Key = Key 
                        
                        j = 0
                        
                        while New_Key in Current_Keys:
                            
                            # Ensure a key is unique - do not overlap data
                            
                            New_Key = Key + str( j )
                            
                            j += 1
                        
                        Compiled_Data[ New_Key ] = Compiled_Data[ Key ]
                        
                        del Compiled_Data[ Key ]
                        
                        New_Keys[ i ] = New_Key
                        
                    else:
                        
                        pass
                    
                #-----------------------------------------------------------------------------------------------------------
                # Add Compiled Data to Pre-Defined Dictionary
                #-----------------------------------------------------------------------------------------------------------
                
                global Count_of_Rows
                
                Count_of_Rows.value = Count_of_Rows.value + len( Compiled_Data )
                
                EQE_Spectra_to_Investigate = EQE_Spectra_to_Investigate | Compiled_Data
        
                EQE_Choice_Widget.value = ''
        
                #-----------------------------------------------------------------------------------------------------------
                # Hide the Data Loading Boxes Now That the Spectra Have Been Added
                #-----------------------------------------------------------------------------------------------------------
    
                Add_Spectrum_Box.layout.display = 'none'
    
                Coda_Box.layout.display = 'none'
                
                Error_Label.layout.display = 'none'
                
                #-----------------------------------------------------------------------------------------------------------
                # Add the Data to the Tab
                #-----------------------------------------------------------------------------------------------------------                
        
                for New_Key in New_Keys:
                
                    EQE_Spectrum_Analysing_Tab_Generator( New_Key )
    else:
        
        Error_Label.layout.display = None
        
        Error_Label.readout = 'Invalid File Name'

# The button executes this function when instructed to by the following code:

Add_Spectrum_Button.on_click( Save_Added_EQE_Spectrum )

# The above function makes use of a data compiling function (defined below) and a pre-existing array that is added to each time the User adds an EQE spectrum.

EQE_Spectra_to_Investigate = {}

def Data_Compiler( Output_DataFrame ):
    
    """When the save data button is pressed, compile the DataFrame into 2D arrays, where each array is stored with one set
    
    of independent and dependent variable data according to the column name."""
    
    Column_Titles = list( Loaded_DataFrame.columns )
    
    #----------------------------------------------------------------------------------------------------------------------
    # Prerequisite Definitions
    #---------------------------------------------------------------------------------------------------------------------- 
        
    Indices_of_Independent_Variables = []
    
    Indices_of_Dependent_Variables = []
    
    Output_Independent_Variable_Data = {}
    
    Output_Dependent_Variable_Data = {}
    
    for Column in Column_Titles:
        
        #------------------------------------------------------------------------------------------------------------------
        # If the column is to be included
        #------------------------------------------------------------------------------------------------------------------
        
        if Include_Checkboxes[ Column ].value == True:   
            
            #--------------------------------------------------------------------------------------------------------------
            # If the data corresponds to an independent variable
            #--------------------------------------------------------------------------------------------------------------

            if Independent_Variable_Checkboxes[ Column ].value == True:

                #----------------------------------------------------------------------------------------------------------
                # Append the column index to the list of independent variable column indices
                #----------------------------------------------------------------------------------------------------------
                
                Indices_of_Independent_Variables.append( Column_Titles.index( Column ) )
                
                #----------------------------------------------------------------------------------------------------------
                # Add the wavelengths to an independent variable data dictionary, convert if necessary
                #----------------------------------------------------------------------------------------------------------
                
                if Type_Widgets[ Column ].value == 'Energy':
                
                    Output_Independent_Variable_Data[ Column ] = Energy_Wavelength_Converter( 
                        
                        array( Loaded_DataFrame[ Column ] ) )
            
                if Type_Widgets[ Column ].value == 'Wavelength':
                
                    Output_Independent_Variable_Data[ Column ] = array( Loaded_DataFrame[ Column ] )
                
            #--------------------------------------------------------------------------------------------------------------
            # If the data corresponds to a dependent variable
            #--------------------------------------------------------------------------------------------------------------

            if Dependent_Variable_Checkboxes[ Column ].value == True:
                
                #----------------------------------------------------------------------------------------------------------
                # Append the column index to the list of dependent variable column indices
                #----------------------------------------------------------------------------------------------------------
                
                Indices_of_Dependent_Variables.append( Column_Titles.index( Column ) )
                
                #----------------------------------------------------------------------------------------------------------
                # Add the EQEs to a dependent variable data dictionary, convert if necessary
                #----------------------------------------------------------------------------------------------------------
                
                if Unit_Widgets[ Column ].value == 'Unitless':
                    
                    Output_Dependent_Variable_Data[ Column ] = array( Loaded_DataFrame[ Column ] )
                    
                if Unit_Widgets[ Column ].value == '%':
                    
                    Output_Dependent_Variable_Data[ Column ] = array( Loaded_DataFrame[ Column ] ) / 100
                    
    #----------------------------------------------------------------------------------------------------------------------
    # Compile Data
    #---------------------------------------------------------------------------------------------------------------------- 
            
    Output_Arrays = {}
    
    for Index in Indices_of_Dependent_Variables:
        
        Preceeding_Independent_Variable_Index = [ Value 
                                                 
                                                 for Value in Indices_of_Independent_Variables if Value < Index ][ -1 ]
                    
        X_Data = Output_Independent_Variable_Data[ Column_Titles[ Preceeding_Independent_Variable_Index ] ]
        
        Y_Data = Output_Dependent_Variable_Data[ Column_Titles[ Index ] ]
        
        if X_Data[ 0 ] > X_Data[ 1 ]:
            
            X_Data = X_Data[ ::-1 ]
            
            Y_Data = Y_Data[ ::-1 ]
        
        Output_Arrays[ Column_Labels[ Column_Titles[ Index ] ].value ] = array( [ X_Data , Y_Data ] )
                        
    return Output_Arrays

# 3.2. Spectral Tailoring

# In this section, the widgets needed to customise irradiance values, irradiance units, superimpose spectra, and more of the like are defined. Firstly, the available spectra are compiled into a series of radio buttons, such that the user can readily switch between them:

Spectrum_Choices = AM_LED_Sheet_Names

Spectrum_Selector = RadioButtons( 
    
    description = 'Spectrum:',
    
    options = AM_LED_Sheet_Names,

    layout={'width': 'max-content'} )

# 3.2.1. Widgets for Superimposing Spectra 

# Following this, widgets are defined for superimposing any number of spectra, starting with checkboxes that "activate" each of the spectra - these are defined and stored in dictionaries using:

Spectrum_Enabled_Boxes = { Spectrum_Type : 
                          
    Checkbox( 
    
        value = False,
    
        description = 'Enable ' + Spectrum_Type,
    
        style = { 'description_width' : 'initial' } )
                         
    for Spectrum_Type in AM_LED_Sheet_Names }

# In addition, input boxes are created for defining the contribution of a given spectrum to the total using:

Spectrum_Ratios_Inputs = { Spectrum_Type : 
                          
    FloatText( 
    
        value = 0,
    
        description = Spectrum_Type + ' Contribution',
        
        style = { 'description_width' : 'initial' },
    
        disabled = True )
                         
    for Spectrum_Type in AM_LED_Sheet_Names }

# Following this, a function is defined for disabling and enabling a given spectrum's ratio input (depending on whether or not the user activates it).

def Spectrum_Input_Disabler( Change ):
    
    """A function for activating/disabling a spectrum's ratio input, depending on whether on what the user chooses. The 
    
    default state is deactivated."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine the spectrum type using the changing widget's description
    #-----------------------------------------------------------------------------------------------------------------------
    
    Spectrum_Type = Change[ 'owner' ].description[ 7: ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # If the widget is activated, enable the input box and give it a value of 1 
    #-----------------------------------------------------------------------------------------------------------------------

    if Change[ 'new' ] == True:
        
        Spectrum_Ratios_Inputs[ Spectrum_Type ].disabled = False
        
        Spectrum_Ratios_Inputs[ Spectrum_Type ].value = 1
        
    #-----------------------------------------------------------------------------------------------------------------------
    # If the widget is deactivated, disable the input box and give it a value of 0
    #-----------------------------------------------------------------------------------------------------------------------
    
    if Change[ 'new' ] == False:
        
        Spectrum_Ratios_Inputs[ Spectrum_Type ].disabled = True
        
        Spectrum_Ratios_Inputs[ Spectrum_Type ].value = 0

# The checkbox widgets are instructed to obey the above function using the following code:

for Spectrum_Type in AM_LED_Sheet_Names:
    
    Spectrum_Enabled_Boxes[ Spectrum_Type ].observe( Spectrum_Input_Disabler , names = 'value' )



# Once the user is happy for their spectra to be superimposed at their desired ratios, they must approve for it to be created using the button defined below.

Approve_Add_Customised_Spectrum_Button = Button(
    
    description = '✓',
    
    button_style = 'success' )

# Alternatively, the user can cancel the creation of their customised spectrum by using the following button:
 
Cancel_Add_Customised_Spectrum_Button = Button(
    
    description = '🗙',
    
    button_style = 'danger' )

# Progress Update:

Progress.value = 85

# However, if the use has not provided ample detail (e.g., only one spectrum that cannot be superimposed with itself), an error message must be revealed. This is done using the following widget:

Spectrum_Customising_Box_Error_Message = Valid(

    value = False,

    readout = 'Invalid',

    style = { 'readout_width' : 'initial' } )

Spectrum_Customising_Box_Error_Message.layout.display = 'None'

# All these widgets for superimposing spectra to create a customised spectrum are compiled into the following widget box:

Spectrum_Customising_Box = VBox( 
    
    [ HBox( [ 
        
        Spectrum_Enabled_Boxes[ Spectrum_Type ] , 
        
        Spectrum_Ratios_Inputs[ Spectrum_Type ] ] )
     
     for Spectrum_Type in AM_LED_Sheet_Names ] +
 
    [ HBox( [ 
        
        Approve_Add_Customised_Spectrum_Button , 
        
        Cancel_Add_Customised_Spectrum_Button ] ),
     
     Spectrum_Customising_Box_Error_Message ] )

# In turn, this box is hidden (until the add customised spectrum process is begun) using

Spectrum_Customising_Box.layout.display = 'None'

# The box should be revealed by pressing the following button:

Add_Customised_Spectrum_Button = Button(

    description = 'Add Customised' )

# The button will do its job using the following function:

def On_Click_Open_Customise_Spectrum_Box( Button ):
    
    """On click, reveral the box containing the widgets for creating customised spectrum."""
    
    Spectrum_Customising_Box.layout.display = None

Add_Customised_Spectrum_Button.on_click( On_Click_Open_Customise_Spectrum_Box )

# The widgets for approving or cancelling the creation of a superimposed spectrum need to be instructed to follow functions to perform their jobs; the latter's function is:

def On_Click_Cancel_Add_Customised_Spectrum( Button ):
    
    """This function cancels the creation of a superimposed spectrum by hiding the customised spectrum creation box, and 
    
    setting all the 'spectrum-enabled' checkboxes to false (which, in turn, sends their values to nought)."""
    
    Spectrum_Customising_Box.layout.display = 'None'
    
    Spectrum_Customising_Box_Error_Message.layout.display = 'None'
    
    for Spectrum_Type in AM_LED_Sheet_Names:
        
        Spectrum_Enabled_Boxes[ Spectrum_Type ].value = False

# The button is now instructed to perform its function using:

Cancel_Add_Customised_Spectrum_Button.on_click( On_Click_Cancel_Add_Customised_Spectrum )

# On the other hand, a function for adding the customised spectrum is defined as follows:

def On_Click_Approve_Add_Customise_Spectrum( Button ):
    
    """On click, add a customised spectrum at the desired ratio, provided that the required criteria are met (i.e., two or
    
    more spectra have been specififed)."""
    
    #----------------------------------------------------------------------------------------------------------------------
    # The enabled spectra are first identified from the user's inputs 
    #---------------------------------------------------------------------------------------------------------------------- 
    
    Desired_Spectra = [ Spectrum_Type for Spectrum_Type in AM_LED_Sheet_Names 
                        
        if Spectrum_Enabled_Boxes[ Spectrum_Type ].value == True ]
    
    #----------------------------------------------------------------------------------------------------------------------
    # If no spectrum has been selected, reveal error message, do nothing
    #----------------------------------------------------------------------------------------------------------------------
    
    if len( Desired_Spectra ) == 0:
        
        Spectrum_Customising_Box_Error_Message.layout.display = None
        
        Spectrum_Customising_Box_Error_Message.readout = 'No spectrum selected'
    
    #----------------------------------------------------------------------------------------------------------------------
    # If only one spectrum has been selected, reveal error message, do nothing
    #----------------------------------------------------------------------------------------------------------------------
        
    if len( Desired_Spectra ) == 1:
        
        Spectrum_Customising_Box_Error_Message.layout.display = None
        
        Spectrum_Customising_Box_Error_Message.readout = 'Select one more'
        
    #----------------------------------------------------------------------------------------------------------------------
    # If two or more spectra have been selected, hide error message (will stay hidden if not revealed), create the spectrum
    #----------------------------------------------------------------------------------------------------------------------
            
    if len( Desired_Spectra ) > 1:
        
        #------------------------------------------------------------------------------------------------------------------
        # Hide error message
        #------------------------------------------------------------------------------------------------------------------
        
        Spectrum_Customising_Box_Error_Message.layout.display = 'None'        
        
        #------------------------------------------------------------------------------------------------------------------        
        # Determine the desired ratios 
        #------------------------------------------------------------------------------------------------------------------
        
        Desired_Ratios = { Spectrum_Type : 
                          
                          Spectrum_Ratios_Inputs[ Spectrum_Type ].value for Spectrum_Type in Desired_Spectra }
        
        #------------------------------------------------------------------------------------------------------------------     
        # Run the external function for superimposing spectra at desired ratios
        #------------------------------------------------------------------------------------------------------------------
        
        Customised_Spectrum_Creator( Desired_Spectra , Desired_Ratios )
        
        #------------------------------------------------------------------------------------------------------------------
        # Hide the box for creating a customised spectrum
        #------------------------------------------------------------------------------------------------------------------
        
        Spectrum_Customising_Box.layout.display = 'None'    

# The button for approving the creation of the superimposed spectra is instructed to obey the above function using the following line of code:

Approve_Add_Customised_Spectrum_Button.on_click( On_Click_Approve_Add_Customise_Spectrum )

# The function relies on another function to superimpose the spectra; namely, the function below:

Additional_Spectra = {}

def Customised_Spectrum_Creator( Desired_Spectra , Desired_Ratios ):
    
    """Create a customised spectrum from a user's specifications."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine the total of the ratios and normalise
    #-----------------------------------------------------------------------------------------------------------------------
        
    Total = sum( list( Desired_Ratios.values() ) )
    
    Normalised_Ratios = { Spectrum_Type : 
                        
                        Desired_Ratios[ Spectrum_Type ] / Total for Spectrum_Type in Desired_Ratios }
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Create a "Combination String" to save the new spectrum with
    #-----------------------------------------------------------------------------------------------------------------------

    Combination_String = ''

    for Key in Desired_Spectra:
    
        Combination_String += Key + ' : ' + str( Desired_Ratios[ Key ] ) + '/' + str( Total ) + ' & '
        
    Combination_String = Combination_String[ :-3 ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Combine the Spectra
    #-----------------------------------------------------------------------------------------------------------------------

    Wavelengths , Combined_Spectrum = Normalised_Photon_Flux_Spectra[ Desired_Spectra[ 0 ] ] 
    
    Combined_Spectrum = Combined_Spectrum * Normalised_Ratios[ Desired_Spectra[ 0 ] ]
    
    for Remaining_Spectrum in Desired_Spectra[ 1: ]:
        
        New_Wavelengths , Normalised_Spectrum = Normalised_Photon_Flux_Spectra[ Remaining_Spectrum ]
        
        Interpolated_Spectrum = Interpolator( New_Wavelengths , Normalised_Spectrum , Wavelengths )
        
        Combined_Spectrum += Interpolated_Spectrum * Normalised_Ratios[ Remaining_Spectrum ]
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Add spectrum to dictionaries
    #-----------------------------------------------------------------------------------------------------------------------
        
    Normalised_Photon_Flux_Spectra[ Combination_String ] = array( [ Wavelengths , Combined_Spectrum ] ) #Spectrum normalised
    
    Photon_Flux_Spectra[ Combination_String ] = array( [ Wavelengths , Combined_Spectrum ] )    
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine P_light and Irradiance Spectra, Then Add Them
    #-----------------------------------------------------------------------------------------------------------------------
    
    Photon_Irradiance_Spectra[ Combination_String ] = array( [ 
        
        Wavelengths , 
        
        10 * Combined_Spectrum * Energy_Wavelength_Converter( Wavelengths ) * e ] )  # Irradiance in units of W/m2/nm
    
    P_lights[ Combination_String ] = simps( 
        
        y = Photon_Irradiance_Spectra[ Combination_String ][ 1 ] , 
                                           
        x = Photon_Irradiance_Spectra[ Combination_String ][ 0 ] ) / 10  # convert to mW/cm2
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Lux Values and Constants of Proportionality, Then Add Them
    #-----------------------------------------------------------------------------------------------------------------------
        
    Arbitrary_Combination_Lux_Values_Calculator( Combination_String )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add spectrum as an option
    #-----------------------------------------------------------------------------------------------------------------------
        
    Spectrum_Selector.options = Spectrum_Selector.options + ( Combination_String , )
    
# The above function creates the superimposed spectrum, adds it as an option to the list of available spectra, stores its photon irradiance spectrum and total spectrum. However, it does not add the lux values and constants of proportionality to the for the already available spectrum. This is done using the function below (which the above function calls):

def Arbitrary_Combination_Lux_Values_Calculator( Spectrum_Label ):
    
    """This function calculates the lux value of the customised spectrum (should be one), before determining the constants
    
    of proportionality that are needed to convert irradiance to lux."""
    
    #----------------------------------------------------------------------------------------------------------------------
    # Load in the spectrum and the wavelengths
    #----------------------------------------------------------------------------------------------------------------------    
    
    Wavelengths , Spectrum = Normalised_Photon_Flux_Spectra[ Spectrum_Label ]
    
    Constants_of_Proportionality_i = {}
    
    Lux_Values_i = {}
    
    #----------------------------------------------------------------------------------------------------------------------
    # Determine the constants of proportionality and the lux values, then append them to their respective dictionaries
    #----------------------------------------------------------------------------------------------------------------------
    
    for Luminous_Efficiency_Type in list( Luminous_Efficiency_Dictionary.keys() ):
    
        Constants_of_Proportionality_i[ Luminous_Efficiency_Type ] = Lux_Constant(
            
            Wavelengths, 
                     
            Photon_Irradiance_Spectra[ Spectrum_Label ][ 1 ] / ( P_lights[ Spectrum_Label ] * 10 ), 
                     
            [ Extrapolated_Luminous_Efficiency( 
                         
                Luminous_Efficiency_Type, 
                         
                Wavelength ) for Wavelength in Wavelengths ] )     
             
        Lux_Values_i[ Luminous_Efficiency_Type ] = Lux_Value( 
            
            Constants_of_Proportionality_i[ Luminous_Efficiency_Type ], 
            
            P_lights[ Spectrum_Label ] * 10  ) # Multiplied by 10 to get W/m2
        
    #----------------------------------------------------------------------------------------------------------------------
    # Output the calculated values
    #----------------------------------------------------------------------------------------------------------------------
    
    Constants_of_Proportionality[ Spectrum_Label ] = Constants_of_Proportionality_i
    
    Lux_Values[ Spectrum_Label ] = Lux_Values_i

# All these widgets for creating a customised are now compiled into one widget box, which can be called into the user interface shortly.

Spectrum_Selector_Box = VBox( [
    
    Spectrum_Selector,
    
    Add_Customised_Spectrum_Button,
    
    Spectrum_Customising_Box ] )

# 3.2.2. Widgets for Customising Light Intensity 

# To specify whether the intensity of the light is given in terms of irradiance or illuminance, the following widget is defined:

Intensity_Type = RadioButtons(
    
    description = 'Intensity Unit',

    options = [ 'Irradiance (W/m2)' , 'Irradiance (mW/cm2)' , 'Illuminance (lx)' ],

    value = 'Irradiance (W/m2)' )

# Similarly, a widget is defined for specifying the type of luminous efficiency spectrum from a few choices, the default it the 2-degree spectrum.

Luminous_Efficiency_Types = list( V_Names.values() )

Luminous_Efficiency_Type_Selector = RadioButtons(
    
    options = Luminous_Efficiency_Types,

    description = 'Luminous Efficiency Type:',

    style = { 'description_width' :'initial' } )

Luminous_Efficiency_Type_Selector.value = 'V 2-deg' 

Luminous_Efficiency_Type_Selector.layout.display = 'none'

# In addition, a widget is defined for specifying the total intensity. The default here is one sun (\sim100\,\mathrm{mW}\cdot\mathrm{cm}^{-2}\approx 116\,\mathrm{k}\,\mathrm{lx}, where one lux is equal to one lumen per square metre, 1\,\mathrm{lx}=1\,\mathrm{lm}\cdot\mathrm{m}^{-2}); but additional options are added to have incrementally sampled customised values (both logarithmically and linear).

Intensity_Sample = RadioButtons(
    
    description = 'Sample Type:',

    options = [ 'One Sun',
               
               'Fixed Value',
               
               'Varied Incrementally (Lin.)',
               
               'Varied Incrementally (Log.)' ] )

# Widgets must now be created such that the user can specify the intensity of light used in their simulations. Firstly, a widget is defined for specifying one customised value:

Intensity_Value = FloatText(
    
    value = 1,
    
    min = 0,

    description = 'Intensity (W/m2)',

    style = { 'description_width' : 'initial' } )

# Additional widgets are defined for incrementally-varied lux values (minimum, maximum, and total number of points; default units are W/m2):

Minimum_Intensity_Value = FloatText(
    
    value = 1,
    
    min = 0,

    description = 'Min. Intensity (W/m2)',

    style = { 'description_width' : 'initial' } )

Maximum_Intensity_Value = FloatText(
    
    value = 100,
    
    min = 0,

    description = 'Max. Intensity (W/m2)',

    style = { 'description_width' : 'initial' } )

Number_of_Points_Value = IntText(
    
    value = 100,
    
    min = 0,

    description = 'Number of Points',

    style = { 'description_width' : 'initial' } )

# These widgets are initially set to be hidden (as one-sun intensity is the default). 

Intensity_Value.layout.display = 'none'

Minimum_Intensity_Value.layout.display = 'none'

Maximum_Intensity_Value.layout.display = 'none'

Number_of_Points_Value.layout.display = 'none'

# Upon changing the unit of intensity, the following widget ensures that the units displayed in the widget labels are correctly updated:

def On_Change_Unit_Changer( Change ):
    
    """If the User changes the unit of the light intensity from irradiance to illuminance, correct the labels of the input 
    
    boxes."""
    
    #---------------------------------------------------------------------------------------------------------------------
    # If the user chooses to use units of W/m2
    #---------------------------------------------------------------------------------------------------------------------    

    if Change[ 'new' ] == 'Irradiance (W/m2)':
        
        Luminous_Efficiency_Type_Selector.layout.display = 'none'
        
        Intensity_Value.description = 'Irradiance (W/m2)'
        
        Minimum_Intensity_Value.description = 'Min. Intensity (W/m2)'
        
        Maximum_Intensity_Value.description = 'Max. Intensity (W/m2)'
                        
        Intensity_Input.description = 'Intensity (W/m2):'
        
        Min_Irradiance.description = 'Minimum Irradiance (W/m2):'

        Max_Irradiance.description = 'Maximum Irradiance (W/m2):'
        
        for Key in list( Intensity_Values.keys() ):
            
            Intensity_Values[ Key ].description = 'Irradiance (W/m2)'
        
    #---------------------------------------------------------------------------------------------------------------------
    # If the user chooses to use units of mW/cm2
    #---------------------------------------------------------------------------------------------------------------------
    
    if Change[ 'new' ] == 'Irradiance (mW/cm2)':
        
        Luminous_Efficiency_Type_Selector.layout.display = 'none'

        Intensity_Value.description = 'Irradiance (mW/cm2)'
        
        Minimum_Intensity_Value.description = 'Min. Intensity (mW/cm2)'
        
        Maximum_Intensity_Value.description = 'Max. Intensity (mW/cm2)'
        
        Intensity_Input.description = 'Intensity (mW/cm2):'
        
        Min_Irradiance.description = 'Minimum Irradiance (mW/cm2):'

        Max_Irradiance.description = 'Maximum Irradiance (mW/cm2):'
        
        for Key in list( Intensity_Values.keys() ):
            
            Intensity_Values[ Key ].description = 'Irradiance (mW/cm2)'        

    #-----------------------------------------------------------------------------------------------------------------------
    # If the user chooses to use units of lux
    #-----------------------------------------------------------------------------------------------------------------------
    
    if Change[ 'new' ] == 'Illuminance (lx)':
        
        Luminous_Efficiency_Type_Selector.layout.display = None
        
        Intensity_Value.description = 'Illuminance (lx)' 
        
        Minimum_Intensity_Value.description = 'Min. Intensity (lx)'
        
        Maximum_Intensity_Value.description = 'Max. Intensity (lx)'
                        
        Intensity_Input.description = 'Intensity (lx):'
        
        Min_Irradiance.description = 'Minimum Illuminance (lx):'

        Max_Irradiance.description = 'Maximum Illuminance (lx):'
        
        for Key in list( Intensity_Values.keys() ):
            
            Intensity_Values[ Key ].description = 'Illuminance (lx)'               
            
# The intensity unit selecting widget is then instructed to observe the above function using:

Intensity_Type.observe( On_Change_Unit_Changer , names = 'value' )

# Furthermore, the appropriate widgets are hidden/revealed according to the user's choice using the following:

def On_Change_Intensity_Customiser( Change ):
    
    """If the User changes the sample type, change the available boxes."""
    
    if Change[ 'new' ] == 'One Sun':
        
        Intensity_Value.layout.display = 'none'

        Minimum_Intensity_Value.layout.display = 'none'

        Maximum_Intensity_Value.layout.display = 'none'

        Number_of_Points_Value.layout.display = 'none'
        
    if Change[ 'new' ] == 'Fixed Value':
        
        Intensity_Value.layout.display = None

        Minimum_Intensity_Value.layout.display = 'none'

        Maximum_Intensity_Value.layout.display = 'none'

        Number_of_Points_Value.layout.display = 'none'

    if Change[ 'new' ] == 'Varied Incrementally (Lin.)' or Change[ 'new' ] == 'Varied Incrementally (Log.)':
        
        Intensity_Value.layout.display = 'none'

        Minimum_Intensity_Value.layout.display = None

        Maximum_Intensity_Value.layout.display = None

        Number_of_Points_Value.layout.display = None

Intensity_Sample.observe( On_Change_Intensity_Customiser , names = 'value' )    

# All these lux customising widgets are then compiled into one box:

Lux_Customisation_Box = VBox( [ 
    
    Spectrum_Selector_Box,
    
    Intensity_Sample,  
        
    HBox( [ 
        
        VBox( [ 
            
            Intensity_Type,
            
            Luminous_Efficiency_Type_Selector ] ),
    
        VBox( [ 
            
            Intensity_Value,
        
            Minimum_Intensity_Value,
            
            Maximum_Intensity_Value,
            
            Number_of_Points_Value ] ) ] ) ] ) 

# Update Progress

Progress.value = 90 

# 3.3. Data Analyser

# In this section, the functions necessary for analysing the data and carrying out the simulations are defined. Using a given input EQE spectrum and irradiance spectrum, the following function will analyser the data to produce values for all figures of merit:

def Data_Analyser( EQE_Spectrum , Light_Spectrum_Label , New_Light_Power, Temperature , Non_Radiative_Loss ):
    
    """For a given EQE spectrum and irradiance spectrum name, calculate the figures of merit. The EQE spectrum is expected 
    
    to be a two dimensional array, with the first row containing wavelength data and the second row containing unitless EQE
    
    values. """
        
    Interpolated_Photon_Flux = Photon_Flux_Interpolator( Light_Spectrum_Label , EQE_Spectrum[ 0 ] , New_Light_Power )
            
    #-----------------------------------------------------------------------------------------------------------------------
    # Short-Circuit Current Density
    #-----------------------------------------------------------------------------------------------------------------------
    
    J_sc = Short_Circuit_Current_Density_Calaculator(
                            
        EQE_Spectrum[ 0 ],
                            
        EQE_Spectrum[ 1 ], 
                            
        Interpolated_Photon_Flux )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Dark Saturation Current
    #-----------------------------------------------------------------------------------------------------------------------
                    
    J_0_rad = Dark_Saturation_Circuit_Current_Density_Calaculator_Rad( 
                        
        EQE_Spectrum[ 0 ], 
                        
        EQE_Spectrum[ 1 ], 
                        
        Planck_Photon_Flux_Wavelength( EQE_Spectrum[ 0 ], Temperature ) )

    J_0 = Dark_Saturation_Circuit_Current_Density_Calaculator( J_0_rad , Non_Radiative_Loss , Temperature )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Radiative Open-Circuit Voltage
    #-----------------------------------------------------------------------------------------------------------------------
                  
    V_oc_rad = Open_Circuit_Voltage( J_sc , J_0_rad , Temperature )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Open-Circuit Voltage
    #-----------------------------------------------------------------------------------------------------------------------
                  
    V_oc = Open_Circuit_Voltage( J_sc, J_0, Temperature )
                    
    #-----------------------------------------------------------------------------------------------------------------------
    # Maximum Power Point Voltage
    #-----------------------------------------------------------------------------------------------------------------------
            
    V_mpp = Maximum_Power_Point_Voltage( V_oc , Temperature )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Maximum Power Point Current
    #-----------------------------------------------------------------------------------------------------------------------
            
    J_mpp = Maximum_Power_Point_Current_Density( V_mpp, J_0, J_sc, Temperature )
                    
    #-----------------------------------------------------------------------------------------------------------------------
    # Maximum Output Power
    #-----------------------------------------------------------------------------------------------------------------------
            
    P_mpp = Maximum_Power_Output( V_mpp, J_mpp )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Fill Factor
    #-----------------------------------------------------------------------------------------------------------------------
            
    FF = Fill_Factor( P_mpp , J_sc , V_oc ) 
                    
    #-----------------------------------------------------------------------------------------------------------------------
    # Power Conversion Efficiency
    #-----------------------------------------------------------------------------------------------------------------------
   
    PCE = Power_Conversion_Efficiency( P_mpp , New_Light_Power )
    
    return { 'J_sc' : J_sc, 
            
            'J_0' : J_0, 
            
            'V_oc' : V_oc, 
            
            'V_oc_rad' : V_oc_rad,
            
            'V_mpp' : V_mpp, 
            
            'J_mpp' : J_mpp, 
            
            'P_mpp' : P_mpp, 
            
            'FF' : FF, 
            
            'PCE' : PCE }

# The following button is defined for analysing the data:

Analyse_Spectra_Button = Button( description = 'Analyse EQE Data' , button_style = 'danger' , disabled = True )

# The colour of the button is instructed to change if the user had added a spectrum:

def On_Change_Analyse_Button_Colour_Changer( Change ):
    
    """Change the colour of the 'Analyse EQE Data' button to give the user an inkling of whether or not the requirements 
    
    have been met."""
    
    if len( EQE_Spectra_to_Investigate ) > 0:
            
        Analyse_Spectra_Button.button_style = 'success'
        
        Analyse_Spectra_Button.disabled = False
            
    else:
            
        Analyse_Spectra_Button.button_style = 'danger'
        
        Analyse_Spectra_Button.disabled = True

EQE_Choice_Widget.observe( On_Change_Analyse_Button_Colour_Changer , names = 'value' )

# The button is instructed to obey the following function:

Analysed_Data = {}

def On_Click_EQE_Spectrum_Analyser( Button ):
    
    """When the user presses the button, analyse the spectra for each of the intensity values."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Firstly, Find the Spectrum Type and the Corresponding Parameters
    #-----------------------------------------------------------------------------------------------------------------------
    
    Spectrum = Spectrum_Selector.value
    
    P_light = P_lights[ Spectrum ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Next, Find the Intensity Sample Type and Unit 
    #-----------------------------------------------------------------------------------------------------------------------
    
    Sample_Type = Intensity_Sample.value
    
    Sample_Unit = Intensity_Type.value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Next, Determine All Intensities to Be Investigated (and Store Them)
    #-----------------------------------------------------------------------------------------------------------------------
    
    if Sample_Unit == 'Irradiance (W/m2)':
        
        Scale_Factor = 1 / 10 
        
    if Sample_Unit == 'Irradiance (mW/cm2)':
        
        Scale_Factor = 1
                
    if Sample_Unit == 'Illuminance (lx)':
        
        V_Type = Luminous_Efficiency_Type_Selector.value
        
        Scale_Factor = 1 / Constants_of_Proportionality[ Spectrum ][ V_Type ] / 10 # Convert to mW/cm2
        
    #-----------------------------------------------------------------------------------------------------------------------
    # If the One Sun Case is To Be Evaluated
    #-----------------------------------------------------------------------------------------------------------------------
        
    if Sample_Type == 'One Sun':
        
        New_Light_Power = P_lights[ 'AM1.5 G']
                
        for Key in list( EQE_Spectra_to_Investigate.keys() ):
                         
            Analysed_Data[ Key ] = { New_Light_Power :
                                        
                Data_Analyser( EQE_Spectra_to_Investigate[ Key ] , Spectrum , New_Light_Power , Temperature_Widget.value , Non_Rad_Loss_Widget.value ) }    
            
    #-----------------------------------------------------------------------------------------------------------------------
    # If a Fixed Value has Been Given
    #-----------------------------------------------------------------------------------------------------------------------
       
    if Sample_Type == 'Fixed Value':
        
        New_Light_Power = Scale_Factor * Intensity_Value.value
                
        for Key in list( EQE_Spectra_to_Investigate.keys() ):
                         
            Analysed_Data[ Key ] = { New_Light_Power :             # Store Data According to Light Power
                                        
                Data_Analyser( EQE_Spectra_to_Investigate[ Key ] , Spectrum , New_Light_Power , Temperature_Widget.value  , Non_Rad_Loss_Widget.value  ) }    
                                                    
    #-----------------------------------------------------------------------------------------------------------------------
    # If Linearly-Spaced Incrementally-Varied Values Have Been Given
    #-----------------------------------------------------------------------------------------------------------------------

    if Sample_Type == 'Varied Incrementally (Lin.)' or Sample_Type == 'Varied Incrementally (Log.)':
        
        min_I = Scale_Factor * Minimum_Intensity_Value.value
        
        max_I = Scale_Factor * Maximum_Intensity_Value.value
        
        N_I = Number_of_Points_Value.value
        
        if Sample_Type == 'Varied Incrementally (Lin.)':
            
            Intensities = linspace( min_I , max_I , int( N_I ) )

        if Sample_Type == 'Varied Incrementally (Log.)':
            
            Intensities = logspace( log10( min_I ) , log10( max_I ) , int( N_I ) )
            
        for Key in list( EQE_Spectra_to_Investigate.keys() ):
                         
            Analysed_Data[ Key ] = { Intensity :
                                        
                Data_Analyser( EQE_Spectra_to_Investigate[ Key ] , Spectrum , Intensity , Temperature_Widget.value , Non_Rad_Loss_Widget.value ) 
                                   
                for Intensity in Intensities }
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Analyse the Data Using the Data Analyser Function
    #-----------------------------------------------------------------------------------------------------------------------
    
Analyse_Spectra_Button.on_click( On_Click_EQE_Spectrum_Analyser )    

# Update Progress

Progress.value = 95 

# 3.4. One-Sun Figures-of-Merit Calculator

# In this section, functions are defined for conducting simulations under one-sun conditions (AM1.5G spectrum at 100\,\mathrm{mW}\,\mathrm{cm}^{-2}). These functions are used to compute the non-radiative open-circuit voltage loss for a given input photovoltaic external quantum efficiency spectrum. This starts with a function that runs the data-analyser function for one-sun conditions: 

def One_Sun_Figures_of_Merit_Calculator( EQE_Spectrum , Temperature ):
    
    """Determine the figures of merit at one sun conditions."""
    
    return Data_Analyser( EQE_Spectrum , 'AM1.5 G' , P_lights[ 'AM1.5 G' ] , Temperature , 0 )

# A function that computes the non-radiative open-circuit voltage loss can then be defined as:

def One_Sun_NR_Loss_Calculator( EQE_Spectrum , Temperature , One_Sun_V_oc ):
    
    """For a given EQE spectrum, calculate the non-radiative loss using the experimental open-circuit voltage at one sun
    
    conditions."""
        
    Interpolated_Photon_Flux = Photon_Flux_Interpolator( 'AM1.5 G' , EQE_Spectrum[ 0 ] , P_lights[ 'AM1.5 G' ] )
            
    #-----------------------------------------------------------------------------------------------------------------------
    # Short-Circuit Current Density
    #-----------------------------------------------------------------------------------------------------------------------
    
    J_sc = Short_Circuit_Current_Density_Calaculator(
                            
        EQE_Spectrum[ 0 ],
                            
        EQE_Spectrum[ 1 ], 
                            
        Interpolated_Photon_Flux )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Radiative Dark Saturation Current
    #-----------------------------------------------------------------------------------------------------------------------
                    
    J_0_rad = Dark_Saturation_Circuit_Current_Density_Calaculator_Rad( 
                        
        EQE_Spectrum[ 0 ], 
                        
        EQE_Spectrum[ 1 ], 
                        
        Planck_Photon_Flux_Wavelength( EQE_Spectrum[ 0 ], Temperature ) )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Radiative Open-Circuit Voltage
    #-----------------------------------------------------------------------------------------------------------------------
                  
    V_oc_rad = k * Temperature / e * log( J_sc / J_0_rad )
    
    return V_oc_rad - One_Sun_V_oc 

# Following this, a function that computes the figures of merit (including the non-radiative open-circuit voltage loss) can be defined as:

def One_Sun_V_oc_Data_Analyser( EQE_Spectrum , Light_Spectrum_Label , New_Light_Power, Temperature , One_Sun_V_oc ):
    
    """Calculate the figures of merit using the V_oc at one sun conditions (AM1.5 G with intensity of 1000 W/m2)."""
    
    Non_Radiative_Loss = One_Sun_NR_Loss_Calculator( EQE_Spectrum , Temperature , One_Sun_V_oc )
    
    if Non_Radiative_Loss < 0:
        
        Non_Radiative_Loss = 0
        
    Figures_of_Merit = Data_Analyser( 
        
        EQE_Spectrum , 
        
        Light_Spectrum_Label , 
        
        New_Light_Power, 
        
        Temperature , 
        
        Non_Radiative_Loss )
    
    Figures_of_Merit[ 'Delta_V_oc_nr' ] = Non_Radiative_Loss
    
    return Figures_of_Merit 

# 4. Supporting Python Tools

# In the following section, additional Python tools that support the simulations are defined, starting with an interpolator function.

# 4.1. Interpolator

# The following function is used to interpolate a data set to estimate what a value would be inbetween data points. It cannot be used to extrapolate:

def Interpolator( x_data , y_data , Desired_x_data ):
    
    """Using a set of x-data and y-data, interpolate to determine what the values would be at the desired x data points.
    
    Do this using NumPy's 'interp' function."""
    
    return interp( Desired_x_data , x_data , y_data )

# In addition, a function is defined for interpolating/extrapolating photon flux measurements. If the measurement lies out of the range of available photon flux data, it is sent to zero:

def Photon_Flux_Interpolator( Light_Spectrum_Type , New_Wavelengths , New_Light_Power ):
    
    """Using a set of x-data and y-data, interpolate to determine what the values would be at the desired x data points."""
    
    #----------------------------------------------------------------------------------------------------------------------
    # Load in the spectral data and total irradiance
    #----------------------------------------------------------------------------------------------------------------------
    
    Wavelengths = Photon_Flux_Spectra[ Light_Spectrum_Type ][ 0 , : ]
    
    Fluxes = Photon_Flux_Spectra[ Light_Spectrum_Type ][ 1 , : ]
    
    P_light = P_lights[ Light_Spectrum_Type ]
    
    #----------------------------------------------------------------------------------------------------------------------
    # Determine bounds on the data
    #----------------------------------------------------------------------------------------------------------------------    
    
    min_wavelength = min( Wavelengths )
    
    max_wavelength = max( Wavelengths )
    
    #----------------------------------------------------------------------------------------------------------------------
    # Define new photon fluxes as the new wavelengths
    #----------------------------------------------------------------------------------------------------------------------    
    
    New_Fluxes = New_Light_Power * array( [ 
        
        Interpolator( Wavelengths, 
                        
                    Fluxes / P_light , 
                        
                    wavelength )
        
        if min_wavelength < wavelength < max_wavelength else 0
        
        for wavelength in New_Wavelengths ] )
    
    return New_Fluxes

# 4.2. Non-Radiative Loss Estimation

# In the following section, non-radiative loss data is imported from an Excel file then parameterised using a few different models. Following this, it is stored in a globally-defined dictionary before being plotted. Storing the parameters of the fitting in a dictionary allows the non-radiative loss to be esimtaed for any optical gap. 

# 4.2.1. Non-Radiative Loss Using Prior Models

# Firstly, the file containing the non-radiative loss data (which is expected to be in different columns of the same file) is specified. 

Non_Radiative_Loss_File_Name = 'Non_Radiative_Loss_Data.xlsx'

Non_Radiative_Loss_Data = read_excel( Non_Radiative_Loss_File_Name )

# The types of non-radiative losses should be given by the column headings of all but the first column:

Non_Radiative_Loss_Types = list( Non_Radiative_Loss_Data.columns )[ 1: ]

# These losses are plotted as a function of the gap using the following code:

# For the smallest optical gaps, these non-radiative losses follow can be modelled by a linear slope of the form \Delta V_oc=mE_\mathrm{g}+c, where m is the gradient of the slope and c is the intercept. A function for modelling this linear graph is defined as:

def Linear_Delta_V_oc( E_g , m , c ) :
    
    """Modelt the open-circuit voltage loss as a straight line."""
    
    return m * E_g + c

# For each of the non-radiative loss types, these gradients and intercepts are now determined using SciPy's 'curve_fit' tool:

Linear_Gradients = {}

Linear_Intercepts = {}

for Loss_Type in Non_Radiative_Loss_Types:
    
    E_gs = list( Non_Radiative_Loss_Data[ Non_Radiative_Loss_Data.columns[ 0 ] ] )[ ::-1 ]   # [ ::-1 ] for increasing E_g 
    
    Delta_V_ocs = list( Non_Radiative_Loss_Data[ Loss_Type ] )[ ::-1 ]
    
    Sampled_E_gs = [ E_gs[ i ] for i in range( len( E_gs ) ) if Delta_V_ocs[ i ] > 0.1 ]
    
    Sampled_Delta_V_ocs = [ Delta_V_oc for Delta_V_oc in Delta_V_ocs if Delta_V_oc > 0.1 ]

    Optimal_Parameters , Covariance = curve_fit(
    
        Linear_Delta_V_oc,
    
        Sampled_E_gs,
    
        Sampled_Delta_V_ocs )
    
    Linear_Gradients[ Loss_Type ] = Optimal_Parameters[ 0 ]
    
    Linear_Intercepts[ Loss_Type ] = Optimal_Parameters[ 1 ]

# The full non-radiative open-circuit voltage loss may be modelled using the following function:

# \Delta V_\mathrm{oc, non-rad}=A\frac{(B-E_\mathrm{g})}{1-\exp(-C[B-E_\mathrm{g}]\right)},

# where A=-m, the gradient of the slope in the linear regime, and B=-\frac{c}{m}. In the limit that the energetic gap is large, Equation ([20](#Delta_V_oc)) reduces to

# \Delta V_\mathrm{oc, non-rad}\approx A(E_\mathrm{g}-B)\exp(-C[E_\mathrm{g}-B]\right)= Ay\exp(-Cy\right),

# where y=E_\mathrm{g}-B. The parameters A and B are determined for each of the loss types below:

As = { Loss_Type : -Linear_Gradients[ Loss_Type ] for Loss_Type in Non_Radiative_Loss_Types }

Bs = { Loss_Type : -Linear_Intercepts[ Loss_Type ]/Linear_Gradients[ Loss_Type ] for Loss_Type in Non_Radiative_Loss_Types }

# The parameters are then used with the following function to determine C. Before this, a function is defined for fitting the data:

def ln_NR_V_oc_Loss_Large_Eg_Approx( y , A , C ):
    
    """Compute the logarithm of the non-radiative open-circuit voltage loss in the large gap approximation."""
    
    return log( A * y ) - C * y

# With the above function defined, the non-radiative open-circuit voltage losses models can be parameterised:

Cs = {}

for Loss_Type in Non_Radiative_Loss_Types:
    
    E_gs = list( Non_Radiative_Loss_Data[ Non_Radiative_Loss_Data.columns[ 0 ] ] )[ ::-1 ]   # [ ::-1 ] for increasing E_g 
    
    Delta_V_ocs = list( Non_Radiative_Loss_Data[ Loss_Type ] )[ ::-1 ]
    
    Sampled_E_gs = [ E_gs[ i ] for i in range( len( E_gs ) ) if Delta_V_ocs[ i ] < 1e-7 ]
    
    Sampled_Delta_V_ocs = [ Delta_V_oc for Delta_V_oc in Delta_V_ocs if Delta_V_oc < 1e-7 ]
    
    A = As[ Loss_Type ]
    
    B = Bs[ Loss_Type ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Apply Curve Fit with Parameter A Fixed to Very Small Bounds (to Prevent it from Varying)
    #-----------------------------------------------------------------------------------------------------------------------

    Small_Bound = 1E-14
    
    Optimal_Parameters , Covariance = curve_fit(
    
        ln_NR_V_oc_Loss_Large_Eg_Approx,
    
        array( Sampled_E_gs ) - B,
    
        log( array( Sampled_Delta_V_ocs ) ),
    
        bounds = ( ( A - Small_Bound , - inf ) , ( A + Small_Bound , inf ) ) )
    
    Cs[ Loss_Type ] = Optimal_Parameters[ -1 ]

# The non-radiative loss parameters are then compiled into a dictionary for each loss type - this can be used to estimate the non-radiative open circuit loss at any bandgap using Equation ([20](#Delta_V_oc)):

Non_Radiative_Loss_Parameters = { Loss_Type : [ As[ Loss_Type ] , Bs[ Loss_Type ] , Cs[ Loss_Type ] ] 
                                 
                                 for Loss_Type in Non_Radiative_Loss_Types }

# The parameters are used by the following function to estimate the non-radiative loss in a given model:

def Non_Radiative_Open_Circuit_Voltage_Loss( E_g , Loss_Type ):
    
    """Determine the open-circuit voltage loss at a given bandgap using the A, B, and C parameter for that loss type."""
    
    A, B, C = Non_Radiative_Loss_Parameters[ Loss_Type ]
    
    return A * ( B - E_g ) / ( 1 - exp( - C * ( B - E_g ) ) )

# A widget for specifying the non-radiative loss is defined below:

Non_Rad_Loss_Widget = FloatText( 
    
    value = 0, 
    
    min = 0, 
    
    description = 'Non-Radiative Loss (V):',

    style = { 'description_width' : 'initial' } )

# 4.2.2. Non-Radiative Loss Using Parabolic Model

# Another realistic model for non-radiative losses based on the data of Ullbrich et al. (Ullbrich, S., et al., _Nature Materials_, __2019__. 18(5): p. 459-464. https://www.nature.com/articles/s41563-019-0324-5; plotted in __Figure 2c__ in the manuscript this computational tool accompanies), 

# \Delta V_\mathrm{oc,non-rad} \approx \begin{cases}
#  0.123 E_\mathrm{g}^2 - 0.64E_\mathrm{g} + 0.927, & \text{ if } E_\mathrm{g}\leq2.60\,\mathrm{eV}.\\
# 0.0945, & \text{ otherwise. }
# \end{cases}

# This equation is encoded using the following Python function:

def Empirical_NR_V_oc_Loss_Ullbrich( E_g ):
    
    """Based on the results of Ullbrich, S., Benduhn, J., Jia, X. et al. Emissive and charge-generating donor–acceptor 
    
    interfaces for organic optoelectronics with low voltage losses. Nat. Mater. 18, 459–464 (2019). 
    
    https://doi.org/10.1038/s41563-019-0324-5"""
    
    if E_g <= 2.601:
    
        return 0.123 * E_g ** 2 - 0.64 * E_g + 0.927 
    
    else: 
        
        return 0.0945

# 4.3. Graph Labels

# In this section, a variety of labels are defined for the output parameters of the simulations. 

Curve_Types = ['J_sc', 'J_0', 'V_oc', 'V_oc_rad', 'V_mpp', 'J_mpp', 'P_mpp', 'FF', 'PCE' ]

# For each of these curve types, a display name is defined (to make the generated graphs more aesthetically appealing):

Display_Names = { 'J_sc' : '$J_\mathrm{sc}$',
                 
               'J_0' : '$J_\mathrm{0}$',
                 
               'V_oc' : '$V_\mathrm{oc}$', 
                 
               'V_oc_rad' : '$V_\mathrm{oc}^\mathrm{rad}$',
                 
               'V_mpp': '$V_\mathrm{mpp}$', 
                 
               'J_mpp': '$J_\mathrm{mpp}$',
                 
               'P_mpp': '$P_\mathrm{mpp}$', 
                 
               'FF' : '$\mathrm{FF}$', 
                 
               'PCE' : '$\mathrm{PCE}$' }

# Moreover, the units of each parameter is defined using:

Curve_Units = { 'J_sc' : '$\mathrm{mA}\cdot\mathrm{cm}^{-2}$',
               
               'J_0' : '$\mathrm{mA}\cdot\mathrm{cm}^{-2}$',
               
               'V_oc' : '$\mathrm{V}$', 
               
               'V_oc_rad' : '$\mathrm{V}$',
               
               'V_mpp': 'V', 
               
               'J_mpp': '$\mathrm{mA}\cdot\mathrm{cm}^{-2}$',
               
               'P_mpp': '$\mathrm{mW}\cdot\mathrm{cm}^{-2}$', 
               
               'FF' : '', 
               
               'PCE' : '' }

# Finally, full labels are defined for each of the parameters (using their display names and units):

Full_Curve_Labels = { Key : Display_Names[ Key ] + ' (' + Curve_Units[ Key ] + ')' for Key in Curve_Types }

Full_Curve_Labels[ 'EQE' ] = 'EQE'

Full_Curve_Labels[ 'Delta_V_oc_nr' ] = '$\\Delta V_\mathrm{oc}^\mathrm{nr} (\mathrm{V})$'

# Following this, lists containing the labels for saving the data (e.g., to Excel files) are specified, starting with the column titles:

Column_Titles = [ 'E_opt', 
                 
                 'J_sc', 
                 
                 'J_0', 
                 
                 'V_oc', 
                 
                 'V_mpp', 
                 
                 'J_mpp', 
                 
                 'P_mpp', 
                 
                 'FF', 
                 
                 'PCE' ]

# Followed by the column units:

Column_Units_Dictionary = { 
    
    'Irradiance (W/m2)' : 'W/m2',
    
    'Irradiance (mW/cm2)' : 'mW/cm2',
    
    'Illuminance (lx)' : 'lx',
    
    'E_opt' : 'eV',
    
    'E_lower' : 'eV',

    'J_sc' : 'mA/cm2',

    'J_0' : 'mA/cm2',

    'V_oc' : 'V',
    
    'V_oc_rad' : 'V',

    'Delta_V_oc_nr' : 'V',

    'V_mpp' : 'V',

    'J_mpp' : 'mA/cm2',

    'P_mpp' : 'mW/cm2',

    'FF' : '',

    'PCE' : '' }

# 4.4. Data Compiler

# In this section, a function is defined for compiling the analysed data (to be used by the simulating tools):

def Analysed_Data_Compiler( Analysed_Data ):
    
    """Compile the data according to lux value."""
    
    EQE_Spectra = list( Analysed_Data.keys() )
    
    Intensity_Values = list( Analysed_Data[ EQE_Spectra[ 0 ] ].keys() )
    
    Compiled_Analysed_Data = { JV_Parameter : 
                              
        { EQE_Spectrum :
            
            [ Analysed_Data[ EQE_Spectrum ][ Intensity ][ JV_Parameter ] for Intensity in Intensity_Values ] 
            
            for EQE_Spectrum in EQE_Spectra }
                              
        for JV_Parameter in Curve_Types }
    
    return Intensity_Values , Compiled_Analysed_Data

# 4.5. Sub-Gap Photovoltaic Quantum Efficiency Simulator 

# In the simulated EQE_PV component of the tool, the photvoltaic quantum efficiency may be modelled in three ways. Firstly, it may be modelled as a step function with some above-gap and below-gap value, it may be modelled as a pseudo-step function (with an exponential tail), or it may be modelled using energetic disorder-dependent exciton absorption (Kaiser, C., et al.  _Nat Commun_ 12, 3988 (__2021__), https://www.nature.com/articles/s41467-021-24202-9; Kay, A., et al., _Adv. Funct. Mater._ __2022__, 32, 2113181, https://onlinelibrary.wiley.com/doi/full/10.1002/adfm.202113181). 

# 4.5.1. Step Function Photovoltaic Quantum Efficiency Simulator 

# The first of these models for the photovoltaic quantum efficiencies is a step function simulator (used to calculate the Shockley-Queisser limit), where all photon of energy greater than some threshold bandgap (E_opt) generate an electron-hole pair at efficiency EQE_{max}, whereas photons of energy less than the threshold bandgap generate at efficiency EQE_\mathrm{min}:

# EQE_PV^\mathrm{SQ}(E)=\begin{cases} 
# EQE_max , \mathrm{\quad if    \quad  } E\geq E_opt,\\
# EQE_\mathrm{min}, \mathrm{\quad otherwise. \quad}
# \end{cases}

# In the Shockley-Queisser limit, the above-gap photovoltaic quantum efficiency EQE_max=1 and the below-gap photovoltaic quantum efficiency EQE_\mathrm{min}=0. These values are left arbitrary here such that the User may customise them. Equation ([23](#SQ_abs)) is encoded using:

def SQ_EQE_Simulator( Energies , Energetic_Gap , Below_Gap_EQE , Above_Gap_EQE ):
    
    """Simulate an EQE sppectrum for a list of energies using values for the above-gap and below-gap EQEs (for the Shockley-
    
    Queisser limit)."""
    
    return array( [ Above_Gap_EQE if E > Energetic_Gap else Below_Gap_EQE for E in Energies ] )

# 4.5.2. Urbach Tail Photovoltaic Quantum Efficiency Simulator 

# Next comes a pseudo-step function, where the sub-gap plateau is replaced with an exponential tail; a more realistic model for absorption described in the work of Urbach (Urbach, F., _Physical Review_, __1953__. 92(5): p. 1324-1324, https://journals.aps.org/pr/abstract/10.1103/PhysRev.92.1324). This exponential tail may be characterised with arbitary Urbach energy E_\mathrm{U}), giving a total photovoltaic quantum efficiency of the form

# EQE_PV^\mathrm{U}(E)=
# EQE_max\begin{cases} 
# 1 , \mathrm{\quad if    \quad  } E\geq E_opt,\\
# \exp(\frac{E-E_opt}{E_\mathrm{U}}\right) , \mathrm{\quad otherwise. \quad}
# \end{cases}

# In organic semiconductors, a reasonable mininum value for the Urbach energy is the thermal energy, E_\mathrm{U}=k_BT (Kaiser, C., et al. A universal Urbach rule for disordered organic semiconductors. _Nat Commun_ 12, 3988 (__2021__). https://doi.org/10.1038/s41467-021-24202-9), where T and k_B were defined in the first section as the temperature and the Boltzmann constant, respectively. In this tool, the temperature is specified using the following widget:

Temperature_Widget = FloatText( value = 293.15,
                               
                               description = 'Temperature, T (K):',
                              
                                style = { 'description_width' : 'initial' } )

# A function for simulating the photovoltaic external quantum efficiency in the sub-gap Urbach tail is defined below:

def E_U_Tail_EQE_Simulator( Energies , Energetic_Gap , Urbach_Energy , Above_Gap_EQE ):
    
    """Simulate an EQE spectrum for a list of energies using values for the above-gap EQE and the Urbach energy."""
    
    return Above_Gap_EQE * array( [ 1 if E > Energetic_Gap else 
                                   
                                   exp( ( E - Energetic_Gap ) / Urbach_Energy ) for E in Energies ] )

# 4.5.3. Organic Semiconductor Photovoltaic Quantum Efficiency Simulator 

# Finally, we have a model describing the absorption of organic semiconductors, in which singlet excitons (SEs) with a disordered density of states (of energetic disorder \sigma_s) and a mean optical gap E_opt have a photovoltaic external quantum efficiency of the form

# EQE_PV^SE(E)=\frac{EQE_max}{2}\{\exp(\frac{E-E_opt+\frac{\sigma_s^2}{2k_BT}}{k_BT}\right)
# erfc(\frac{E-E_opt+\frac{\sigma_s^2}{k_BT}}{\sigma_s\sqrt{2}}\right)+
# erf(\frac{E_opt}{\sigma_s\sqrt{2}}\right)+erf(\frac{E-E_opt}{\sigma_s\sqrt{2}}\right)\right\}
# \approx\frac{EQE_max}{2}\{\exp(\frac{E-E_opt+\frac{\sigma_s^2}{2k_BT}}{k_BT}\right)
# erfc(\frac{E-E_opt+\frac{\sigma_s^2}{k_BT}}{\sigma_s\sqrt{2}}\right)+
# erfc(\frac{E_opt-E}{\sigma_s\sqrt{2}}\right)\right\},

# where erf and erfc denote the error function and complementary error function, respectively. This equation is defined using the following function:

from scipy.special import erf, erfc

def SE_EQE_Simulator( Energies , Energetic_Gap , Energetic_Disorder , Above_Gap_EQE , Temperature ):
    
    """Simulate an EQE spectrum for a list of energies using values for the above-gap EQE, the energetic gap, the energetic
    
    disorder, and the temperature."""
    
    kT = k * Temperature / e
    
    if Energetic_Disorder != 0:
        
        Exponential_Term = exp( ( Energies - Energetic_Gap + Energetic_Disorder ** 2 / 2 / kT  ) / kT )

        Erfc_Term = erfc( ( Energies - Energetic_Gap + Energetic_Disorder ** 2 / kT  ) / Energetic_Disorder / 2 ** 0.5 )
        
        Erf_Term_1 = erf( ( Energetic_Gap ) / Energetic_Disorder / 2 ** 0.5 )
        
        Erf_Term_2 = erf( ( Energies - Energetic_Gap ) / Energetic_Disorder / 2 ** 0.5 ) 
        
        return Above_Gap_EQE / 2 * ( Exponential_Term * Erfc_Term + Erf_Term_1 + Erf_Term_2 )
        
    else:
        
        return E_U_Tail_EQE_Simulator( Energies , Energetic_Gap , kT , Above_Gap_EQE )
#  4.6 Data Saving and Copying Tools

# Throughout this section, the tools needed to compile and save the data generated by the simulations are defined, starting with the tools for the optical gap-dependent simulations in [Section 4.6.1](#E_opt_Data_Saving) and intensity-dependent simulations in [Section 4.6.2](#I_Data_Saving). Firstly, some functions are defined for compiling the figures of merit in each case:

def Figures_of_Merit_Compiler( Energetic_Gaps , Figures_of_Merit_Dictionary ):
    
    """For energetic gap-dependent simulations, compile the figures of merit."""
    
    Compiled_Figures_of_Merit = { 'E_opt' : list( Energetic_Gaps ) }
    
    for Curve_Type in Curve_Types:
        
        Compiled_Figures_of_Merit[ Curve_Type ] = [
            
            Figures_of_Merit_Dictionary[ Gap ][ Curve_Type ] for Gap in Energetic_Gaps ]
    
    return Compiled_Figures_of_Merit

def Figures_of_Merit_Compiler_vs_Lux( Intensities , Figures_of_Merit_Dictionary ):
    
    """For intensity-dependent simulations, compile the figures of merit."""
    
    Compiled_Figures_of_Merit = { 'Intensities' : list( Intensities ) }
    
    for Curve_Type in Curve_Types:
        
        Compiled_Figures_of_Merit[ Curve_Type ] = [
            
            Figures_of_Merit_Dictionary[ Intensity ][ Curve_Type ] for Intensity in Intensities ]
    
    return Compiled_Figures_of_Merit

# 4.6.1. Saving Data for Optical Gap-Dependent Simulations

# To save the data, a Pandas data frame must be created because this can be saved directly to an Excel file. This file path must first be generated; the following imports are therefore made:

from os import mkdir 
from pandas import DataFrame

# The following function is now defined for saving the figure-of-merit versus optical gap data to an excel file:

def On_Click_Figure_of_Merit_vs_E_opt_Saver( Button ):
    
    """On click, save the compiled figures of merit dictionary as a text document in a Folder."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Output Directory
    #-----------------------------------------------------------------------------------------------------------------------
    
    Output_Directory_Path = path.join( Current_Working_Directory , 'Output_Data' )
    
    if 'Output_Data' not in listdir( Current_Working_Directory ):
        
        mkdir( Output_Directory_Path )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Output Directory (to be converted to a dataframe) and Export Data
    #-----------------------------------------------------------------------------------------------------------------------
    
    Loss_Types_to_Export = list( Compiled_Figures_of_Merit.keys() )

    for Loss_Type in Loss_Types_to_Export:
        
        File_Name = 'Figures_of_Merit_vs_E_opt_' + Loss_Type + '.txt'
        
        File_Path = path.join( Output_Directory_Path , File_Name )
        
        Output_DataFrame = DataFrame( Compiled_Figures_of_Merit[ Loss_Type ] )
    
        #-------------------------------------------------------------------------------------------------------------------
        # Add Column Titles and Units 
        #-------------------------------------------------------------------------------------------------------------------
    
        Title_DataFrame = DataFrame( data = None , 
                                    
                                    columns = Output_DataFrame.columns )
        
        Title_DataFrame.loc[ 0 ] = Output_DataFrame.columns 
        
        Title_DataFrame.loc[ 1 ] = [ Column_Units_Dictionary[ Key ] for Key in Output_DataFrame.columns ] 
        
        with open( File_Path , 'w' ) as Output_File:
            
            savetxt( Output_File , 
                 
                 Title_DataFrame.values , 
                
                 fmt = '%s',
                                 
                 delimiter = '\t', 
                 
                 newline = '\n' )
            
        Output_File.close()
        
        #-------------------------------------------------------------------------------------------------------------------
        # Save Data
        #-------------------------------------------------------------------------------------------------------------------
    
        with open( File_Path , 'a' ) as Output_File:
            
            savetxt( Output_File , 
                 
                 Output_DataFrame.values , 
                                 
                 delimiter = '\t', 
                 
                 newline = '\n' )
            
        Output_File.close()
        
# 4.6.2. Saving Data for Intensity-Dependent Simulations

# In a similar way to the above function, the following function will save the figures-of-merit for intensity-dependent simulations:

def On_Click_Figure_of_Merit_vs_Intensity_Saver( Button ):
    
    """On click, save the compiled figures of merit dictionary as a text document in a Folder."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Output Directory
    #-----------------------------------------------------------------------------------------------------------------------
    
    Output_Directory_Path = path.join( Current_Working_Directory , 'Output_Data' )
    
    if 'Output_Data' not in listdir( Current_Working_Directory ):
        
        mkdir( Output_Directory_Path )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Output Directory and Export Data
    #-----------------------------------------------------------------------------------------------------------------------
    
    Loss_Types_to_Export = list( Compiled_Figures_of_Merit_vs_Intensity.keys() )

    for Loss_Type in Loss_Types_to_Export:
        
        File_Name = 'Figures_of_Merit_vs_Intensity_' + Loss_Type + '.txt'
        
        File_Path = path.join( Output_Directory_Path , File_Name )
        
        Output_DataFrame = DataFrame( Compiled_Figures_of_Merit_vs_Intensity[ Loss_Type ] )
    
        #-------------------------------------------------------------------------------------------------------------------
        # Add Column Titles and Units 
        #-------------------------------------------------------------------------------------------------------------------
    
        Title_DataFrame = DataFrame( data = None , 
                                    
                                    columns = Output_DataFrame.columns )
        
        Title_DataFrame.loc[ 0 ] = Output_DataFrame.columns 
        
        Title_DataFrame.loc[ 1 ] = [ Column_Units_Dictionary[ Intensity_Type.value ] ] + [ 
            
            Column_Units_Dictionary[ Key ] for Key in list( Output_DataFrame.columns )[ 1: ] ] 
        
        savetxt( File_Path , 
                 
                 Title_DataFrame.values , 
                
                 fmt = '%s',
                                 
                 delimiter = '\t', 
                 
                 newline = '\n' )
        
        #-------------------------------------------------------------------------------------------------------------------
        # Save Data
        #-------------------------------------------------------------------------------------------------------------------
    
        with open( File_Path , 'a' ) as Output_File:
            
            savetxt( Output_File , 
                 
                 Output_DataFrame.values , 
                                 
                 delimiter = '\t', 
                 
                 newline = '\n' )
            
        Output_File.close()

# 4.6.3. Copying Data for Optical Gap-Dependent Simulations

# Following this, functions are defined for optical gap-dependent copying data (to the User's clipboard), allowing the data to be readily pasted into another program (like, e.g., Excel). These functions rely on the following function for copying Python arrays to the clipboard:

import win32clipboard as clipboard

def toClipboard( Array ):
    
    """Copies an array into a string format acceptable by Excel. Columns separated by \t, rows separated by \n. From
    
    https://stackoverflow.com/questions/22488566/how-to-paste-a-numpy-array-to-excel"""
    
    # Create string from array
    
    Line_Strings = []
    
    for line in Array:
    
        Line_Strings.append( "\t".join( line.astype( str ) ).replace( "\n" , "" ) )
    
    Array_String = "\r\n".join( Line_Strings )

    # Put string into clipboard (open, clear, set, close)
    
    clipboard.OpenClipboard()
    
    clipboard.EmptyClipboard()
    
    clipboard.SetClipboardText( Array_String )
    
    clipboard.CloseClipboard()

# The following function can then be used to copy optical gap-dependent data to the clipboard (for the simulated EQE_PV case). It identifies how many loss types the user has enabled, before creating an output array containing the headers, the units, the loss types, and the data (including the optical gap data): 

def On_Click_E_opt_Dep_Data_to_Clipboard_Copier( Button ):
    
    """Copy E_opt-dependent data to clipboard."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Prepare Output Array
    #-----------------------------------------------------------------------------------------------------------------------    
    
    Loss_Types = list( Compiled_Figures_of_Merit.keys() ) 
    
    Number_of_Entries = len( Compiled_Figures_of_Merit[ Loss_Types[ 0 ] ].keys() )
    
    Array_Width = Number_of_Entries * len( Loss_Types )
    
    Array_Length = len( Compiled_Figures_of_Merit[ Loss_Types[ 0 ] ][ 'E_opt' ] )
    
    Output_Array = zeros( [ Array_Length , Array_Width ] )
    
    Header_Array = []
    
    Title_Array = []
    
    Unit_Array = []
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Populate Output Array
    #-----------------------------------------------------------------------------------------------------------------------        
    
    for i in range( len( Loss_Types ) ): 
        
        Data_Set_i = Compiled_Figures_of_Merit[ Loss_Types[ i ] ]
        
        Keys_i = list( Data_Set_i.keys() )
        
        for j in range( Number_of_Entries ):
            
            Header_Array.append( Loss_Types[ i ] )
            
            Key_j = Keys_i[ j ]
            
            Title_Array.append( Key_j )
            
            Unit_Array.append( Column_Units_Dictionary[ Key_j ] )
            
            Output_Array[ : , Number_of_Entries * i + j ] = Data_Set_i[ Key_j ]
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Copy Output Array to Clipboard
    #-----------------------------------------------------------------------------------------------------------------------    
        
    New_Output_Array = concatenate( ( 
        
        array( [ Title_Array , Unit_Array , Header_Array ] ),
        
        Output_Array ) )
    
    toClipboard( New_Output_Array )

# On the other hand, the function that comes next copies optical gap-dependent data to the clipboard in the case of experimental  EQE_PV. It identifies how many loss types the user has enabled, before creating an output array containing the headers, the units, the loss types, and the data (including the optical gap data): 

def On_Click_EQE_Versus_E_Lower_to_Clipboard_Copier( Button ):
    
    """Copy E_opt-dependent data to clipboard."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Prepare Output Array
    #-----------------------------------------------------------------------------------------------------------------------        
        
    Key = Button.style._view_name
        
    Limit_Dep_Figures_of_Merit = Limit_Dep_Figures_of_Merits[ Key ]

    Column_Titles = list( Limit_Dep_Figures_of_Merit.keys() )
    
    Array_Width = len( Column_Titles )
    
    Array_Length = len( Limit_Dep_Figures_of_Merit[ Column_Titles[ 0 ] ] )
    
    Output_Array = zeros( [ Array_Length , Array_Width ] )
    
    Header_Array = []
    
    Title_Array = []
    
    Unit_Array = []
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Populate Output Array
    #-----------------------------------------------------------------------------------------------------------------------            
    
    for j in range( Array_Width ):
        
        Column_Title = Column_Titles[ j ]
            
        Header_Array.append( Key )
            
        Title_Array.append( Column_Title )
        
        if Column_Titles[ j ] == 'Intensity':
            
            Unit_Array.append( Intensity_Type.value )
            
        else:
            
            Unit_Array.append( Column_Units_Dictionary[ Column_Title ] )
            
        Output_Array[ : , j ] = Limit_Dep_Figures_of_Merit[ Column_Title ]
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Copy Output Array to Clipboard
    #-----------------------------------------------------------------------------------------------------------------------        
            
    New_Output_Array = concatenate( ( 
        
        array( [ Title_Array , Unit_Array , Header_Array ] ),
        
        Output_Array ) )
    
    toClipboard( New_Output_Array )
    
# 4.6.4. Copying Data for Intensity-Dependent Simulations

# In this sub-section, functions are defined for copying intensity-dependent data (to the User's clipboard), allowing the data to be readily pasted into another program (like, e.g., Excel). These functions rely on the following function for copying Python arrays to the clipboard:

def On_Click_Intensity_Dep_Data_to_Clipboard_Copier( Button ):
    
    """Copy Intensity-dependent data to clipboard."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Prepare Output Array
    #-----------------------------------------------------------------------------------------------------------------------            
    
    Loss_Types = list( Compiled_Figures_of_Merit_vs_Intensity.keys() ) 
    
    Number_of_Entries = len( Compiled_Figures_of_Merit_vs_Intensity[ Loss_Types[ 0 ] ].keys() )
    
    Array_Width = Number_of_Entries * len( Loss_Types )
    
    Array_Length = len( Compiled_Figures_of_Merit_vs_Intensity[ Loss_Types[ 0 ] ][ 'Intensities' ] )
    
    Output_Array = zeros( [ Array_Length , Array_Width ] )
    
    Header_Array = []
    
    Title_Array = []
    
    Unit_Array = []
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Populate Output Array
    #-----------------------------------------------------------------------------------------------------------------------            
    
    for i in range( len( Loss_Types ) ): 
        
        Data_Set_i = Compiled_Figures_of_Merit_vs_Intensity[ Loss_Types[ i ] ]
        
        Keys_i = list( Data_Set_i.keys() )
        
        for j in range( Number_of_Entries ):
            
            Header_Array.append( Loss_Types[ i ] )
            
            Key_j = Keys_i[ j ]
                
            Title_Array.append( Key_j )
            
            if Key_j == 'Intensities':
            
                Unit_Array.append( Column_Units_Dictionary[ Intensity_Type.value ] )
                
            else:
                
                Unit_Array.append( Column_Units_Dictionary[ Key_j ] )
            
            Output_Array[ : , Number_of_Entries * i + j ] = Data_Set_i[ Key_j ]
            
    #-----------------------------------------------------------------------------------------------------------------------
    # Copy Output Array to Clipboard
    #-----------------------------------------------------------------------------------------------------------------------                    
    
    New_Output_Array = concatenate( ( 
        
        array( [ Title_Array , Unit_Array , Header_Array ] ),
        
        Output_Array ) )
    
    toClipboard( New_Output_Array )  

# Where the above function copied intensity-dependent data in the simulated EQE_PV case, the below function copies data in the experimental EQE_PV case:

def On_Click_EQE_Versus_Intensity_to_Clipboard_Copier( Button ):
    
    """Copy E_opt-dependent data to clipboard."""
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Prepare Output Array
    #-----------------------------------------------------------------------------------------------------------------------                
    
    Key = Button.style._view_name
    
    Intensity_Dep_Figures_of_Merit = Intensity_Dep_Figures_of_Merits[ Key ]

    Column_Titles = list( Intensity_Dep_Figures_of_Merit.keys() )
    
    Array_Width = len( Column_Titles )
    
    Array_Length = len( Intensity_Dep_Figures_of_Merit[ Column_Titles[ 0 ] ] )
    
    Output_Array = zeros( [ Array_Length , Array_Width ] )
    
    Header_Array = []
    
    Title_Array = []
    
    Unit_Array = []
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Populate Output Array
    #-----------------------------------------------------------------------------------------------------------------------            
    
    for j in range( Array_Width ):
        
        Column_Title = Column_Titles[ j ]
            
        Header_Array.append( Key )
            
        Title_Array.append( Column_Title )
        
        if Column_Titles[ j ] == 'Intensity':
            
            Unit_Array.append( Intensity_Type.value )
            
        else:
            
            Unit_Array.append( Column_Units_Dictionary[ Column_Title ] )
            
        Output_Array[ : , j ] = Intensity_Dep_Figures_of_Merit[ Column_Title ]
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Copy Output Array to Clipboard
    #-----------------------------------------------------------------------------------------------------------------------                    
        
    New_Output_Array = concatenate( ( 
        
        array( [ Title_Array , Unit_Array , Header_Array ] ),
        
        Output_Array ) )
    
    toClipboard( New_Output_Array )        

# 4.7. Control Functions

# The functions defined in this section are, by far, the most involved functions of the lot. They control the simulations, update the interface, and store the data. This section starts with the functions that control the simulated EQE_PV case in [Section 4.7.1.](#ControlFunctionSimEQE), which is followed by the functions that control the experimental EQE_PV case in [Section 4.7.2.](#ControlFunctionExpEQE). 

# 4.7.1. Controlling in the Simulated EQE Case

# The following two functions conduct simulate the figures of merit in the simulated EQE_PV case, before storing the data, generating the graphs, and updating the interface. These functions tie together all the other functions of the script - they are therefore the biggest. The first function controls simulations versus the optical gap:

def On_Click_Limit_vs_E_opt_Computer( Button ):
    
    """Compute the limits versus optical gap on click of the correct button."""

    #----------------------------------------------------------------------------------------------------------------------
    # Load in Photon Flux Spectrum
    #----------------------------------------------------------------------------------------------------------------------
    
    Spectrum_Type = Spectrum_Selector.value
    
    Energies = linspace( Min_Sim_Energy.value , Max_Sim_Energy.value , N_Sim_Energies.value )[ ::-1 ]
    
    Wavelengths = Energy_Wavelength_Converter( Energies )
    
    #----------------------------------------------------------------------------------------------------------------------
    # Determine Bandgaps
    #----------------------------------------------------------------------------------------------------------------------
    
    Calculation_Types = [ Type for Type in PCE_Calculation_Types if PCE_Calculation_Select_Widgets[ Type ].value == True ] 
    
    Energetic_Gaps = linspace( Min_Energetic_Gap.value , Max_Energetic_Gap.value , N_Energetic_Gaps.value )
    
    global Compiled_Figures_of_Merit
    
    Compiled_Figures_of_Merit = {}
    
    Peak_Performance_Labels = []
    
    for Type in Calculation_Types:
                
        if EQE_Simulation_Type.value == 'Step Function':
            
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Step Function Like
            #--------------------------------------------------------------------------------------------------------------
                
            EQE_Spectra = { Energetic_Gap : 
                           
                array( [ Wavelengths,
                    
                    SQ_EQE_Simulator( Energies, 
                                                             
                        Energetic_Gap, 
                                                             
                        Below_Gap_Value.value,
                                                            
                        Above_Gap_Value.value ) ] )
                           
                for Energetic_Gap in Energetic_Gaps }
        
        if EQE_Simulation_Type.value == 'Urbach Tail (Inorganic/Perovskite)':
        
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Urbach Tail Discontinuous Transition
            #--------------------------------------------------------------------------------------------------------------
                
            if Use_Thermal_Energy_Checkbox.value == True:
                
                    E_U = Temperature_Widget.value * k / e
                    
            else:
                
                E_U = Urbach_Energy_Value.value
                
            EQE_Spectra = { Energetic_Gap : 
                           
                array( [ Wavelengths,
                    
                    E_U_Tail_EQE_Simulator( Energies, 
                                                             
                        Energetic_Gap, 
                                                             
                        E_U,
                                                            
                        Above_Gap_Value.value ) ] )
                           
                for Energetic_Gap in Energetic_Gaps }
            
        if EQE_Simulation_Type.value == 'Exciton Absorption (Organic)':
        
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Excitonic Absorption
            #--------------------------------------------------------------------------------------------------------------
                
            EQE_Spectra = { Energetic_Gap : 
                           
                array( [ Wavelengths,
                    
                    SE_EQE_Simulator( Energies , 
                                     
                        Energetic_Gap, 
                                     
                        Energetic_Disorder_Value.value, 
                                     
                        Above_Gap_Value.value , 
                                     
                        Temperature_Widget.value ) ] )
                           
                for Energetic_Gap in Energetic_Gaps }
        
        #--------------------------------------------------------------------------------------------------------------
        # Determine Non-Radiative Loss
        #--------------------------------------------------------------------------------------------------------------
            
        if Type == 'Radiative Limit':
                
            V_oc_Losses = { Energetic_Gap : 0 for Energetic_Gap in Energetic_Gaps }
            
        else:
        
            if Type == 'Empirical - Ullbrich et al.':
            
                V_oc_Losses = { Energetic_Gap : 
                           
                    Empirical_NR_V_oc_Loss_Ullbrich( Energetic_Gap ) for Energetic_Gap in Energetic_Gaps  }
                
            else:
                
                V_oc_Losses = { Energetic_Gap : 
                               
                    Non_Radiative_Open_Circuit_Voltage_Loss( Energetic_Gap , Type )
                               
                    for Energetic_Gap in Energetic_Gaps }
                
        #--------------------------------------------------------------------------------------------------------------
        # Determine New Light Power
        #--------------------------------------------------------------------------------------------------------------

        if Intensity_Type.value == 'Irradiance (W/m2)':
            
            Scale_Factor = 1
            
        if Intensity_Type.value == 'Irradiance (mW/cm2)':
            
            Scale_Factor = 0.1
            
        if Intensity_Type.value == 'Illuminance (lx)':
            
            Scale_Factor = P_lights[ Spectrum_Type ] / Lux_Values[ Spectrum_Type ][ 'V 2-deg' ]
          
        New_Light_Power = Scale_Factor * Intensity_Input.value
        
        #--------------------------------------------------------------------------------------------------------------
        # Determine Figures of Merit
        #--------------------------------------------------------------------------------------------------------------
        
        Figures_of_Merit = { Energetic_Gap : 
                                
            Data_Analyser( EQE_Spectra[ Energetic_Gap ],
                              
                Spectrum_Type, 
                              
                New_Light_Power,
                          
                Temperature_Widget.value, 
                             
                V_oc_Losses[ Energetic_Gap ] )
                                
            for Energetic_Gap in Energetic_Gaps }
                    
        Compiled_Figures_of_Merit_Local = Figures_of_Merit_Compiler( Energetic_Gaps , Figures_of_Merit )
            
        Compiled_Figures_of_Merit[ Type ] = Compiled_Figures_of_Merit_Local
            
        #--------------------------------------------------------------------------------------------------------------
        # Determine Maximum PCE
        #--------------------------------------------------------------------------------------------------------------
    
        PCEs = Compiled_Figures_of_Merit_Local[ 'PCE' ]
            
        Optimal_Index = where( PCEs == max( PCEs ) )[ 0 ][ 0 ] 
            
        Optimal_Gap = Energetic_Gaps[ Optimal_Index ] 
            
        Loss_Type_Label = Label( value = Type + ':' )
        
        Optimal_Gap_Label = Label( value = 'Optimal Gap = ' + str( Optimal_Gap ) + ' eV')
        
        Max_PCE_Label = Label( value = 'Max PCE = ' + str( PCEs[ Optimal_Index ] ) )
        
        Peak_Performance_Labels.append( 
            
            VBox( [ 
        
                Loss_Type_Label,
            
                Optimal_Gap_Label,
            
                Max_PCE_Label ] ) )
        
    #--------------------------------------------------------------------------------------------------------------
    # Generate Graphs
    #--------------------------------------------------------------------------------------------------------------
    
    global Output_Graphs
    
    Output_Graphs = { Curve_Type : Output() for Curve_Type in Curve_Types }
    
    for Curve_Type in Curve_Types:
                
        if Curve_Units[ Curve_Type ] != '':
                    
            Unit = ' (' + Curve_Units[ Curve_Type ] + ')'
                    
        else:
                    
            Unit = ''
                
        with Output_Graphs[ Curve_Type ]:
                    
            plt.figure( dpi = 100 )
        
            for Type in Calculation_Types:
                    
                plt.plot( Energetic_Gaps , Compiled_Figures_of_Merit[ Type ][ Curve_Type ] , label = Type )
                                        
            plt.ylabel( Display_Names[ Curve_Type ] + Unit  )
            
            plt.xlabel( '$E_\mathrm{g}$ (eV)')
                    
            if Curve_Type == 'J_0':
                        
                plt.yscale( 'log' )
            
            plt.legend()
                    
            plt.show()
                    
        Output_Graphs[ Curve_Type ]
                    
        Output_Graphs[ Curve_Type ].layout.display = 'none'

    #--------------------------------------------------------------------------------------------------------------
    # Output Graphs and Enable Save Data Button
    #--------------------------------------------------------------------------------------------------------------
    
    Save_E_opt_Dep_Data_Button.disabled = False
    
    Copy_Limit_vs_E_opt_Data_to_Clipboard_Button.disabled = False
        
    #Output_Graphs[ Curve_Types[ 0 ] ].layout.display = None
            
    Curve_Control = RadioButtons( 
                
        options = Curve_Types,
            
        value = Curve_Types[ 0 ] )
            
    Compiled_Box = VBox( 
        
        [ Curve_Control ] + 
        
        [ Output_Graphs[ Curve_Type ] for Curve_Type in Curve_Types ] +
    
        [ Label_Box for Label_Box in Peak_Performance_Labels ] )
            
    Curve_Control.observe( On_Change_Output_Graph_Revealer, names = 'value' )
            
    Curve_Control.value = 'PCE'
                
    Varied_E_opt_Box.children = Varied_E_opt_Box.children[ :-1 ]
            
    Varied_E_opt_Box.children = Varied_E_opt_Box.children + ( Compiled_Box , )

# This first function relies on one further external function to change which graph is currently revealed on the interface:

def On_Change_Output_Graph_Revealer( Change ):
    
    """Change which graph is revealed in the simulation pane.""" 
    
    for Curve_Type in Curve_Types:
        
        Output_Graphs[ Curve_Type ].layout.display = 'none'
        
    Output_Graphs[ Change[ 'new' ] ].layout.display = None

# The second function, on the other hand, controls simulations versus intensity:

def On_Click_Limit_vs_Intensity_Computer( Button ):
    
    """Compute the limits versus intensity on click of the correct button."""

    #----------------------------------------------------------------------------------------------------------------------
    # Load in Photon Flux Spectrum
    #----------------------------------------------------------------------------------------------------------------------
    
    Spectrum_Type = Spectrum_Selector.value
    
    Energies = logspace( log10( Min_Sim_Energy.value ) , log10( Max_Sim_Energy.value ) , N_Sim_Energies.value )[ ::-1 ]
    
    Wavelengths = Energy_Wavelength_Converter( Energies )
    
    #----------------------------------------------------------------------------------------------------------------------
    # Find Optimal Gap
    #----------------------------------------------------------------------------------------------------------------------
    
    Calculation_Types = [ Type for Type in PCE_Calculation_Types if PCE_Calculation_Select_Widgets[ Type ].value == True ] 
    
    if Determine_Optimal_Gap_Checkbox.value == False:
    
        Energetic_Gaps = { Type : E_opt_Value.value
                          
                          for Type in Calculation_Types }
        
    else:
        
        Decimal_Places = int( E_opt_Value.value )
        
        Energetic_Gaps = {}
        
        for Type in Calculation_Types:
            
            Best_E_opt = 1.5
        
            for i in range( Decimal_Places ):
            
                Scale_Factor = 1 / 10 ** i      # Scale factor is 1 is i = 0, 1/10 if i = 1 , etc.
            
                Sample_E_opts = linspace( - 1 , 1 , 21 ) * Scale_Factor + Best_E_opt 
                
                Sample_E_opts = array( [ E_opt for E_opt in Sample_E_opts if E_opt > 0 ] ) # Remove any negative gaps
                
                if EQE_Simulation_Type.value == 'Step Function':
            
                    #------------------------------------------------------------------------------------------------------
                    # Simulate Spectra - Step Function Like
                    #------------------------------------------------------------------------------------------------------
                
                    EQE_Spectra = { Energetic_Gap : 
                           
                        array( [ Wavelengths,
                    
                            SQ_EQE_Simulator( Energies, 
                                                             
                                Energetic_Gap, 
                                                             
                                Below_Gap_Value.value,
                                                            
                                Above_Gap_Value.value ) ] )
                           
                        for Energetic_Gap in Sample_E_opts }
        
                if EQE_Simulation_Type.value == 'Urbach Tail (Inorganic/Perovskite)':
        
                    #------------------------------------------------------------------------------------------------------
                    # Simulate Spectra - Urbach Tail Discontinuous Transition
                    #------------------------------------------------------------------------------------------------------
                
                    if Use_Thermal_Energy_Checkbox.value == True:
                
                        E_U = Temperature_Widget.value * k / e
                    
                    else:
                
                        E_U = Urbach_Energy_Value.value
                
                    EQE_Spectra = { Energetic_Gap : 
                           
                        array( [ Wavelengths,
                    
                            E_U_Tail_EQE_Simulator( Energies, 
                                                             
                                Energetic_Gap, 
                                                             
                                E_U,
                                                            
                                Above_Gap_Value.value ) ] )
                           
                        for Energetic_Gap in Sample_E_opts }
            
                if EQE_Simulation_Type.value == 'Exciton Absorption (Organic)':
        
                    #------------------------------------------------------------------------------------------------------
                    # Simulate Spectra - Excitonic Absorption
                    #------------------------------------------------------------------------------------------------------
                
                    EQE_Spectra = { Energetic_Gap : 
                           
                        array( [ Wavelengths,
                    
                            SE_EQE_Simulator( Energies , 
                                     
                            Energetic_Gap, 
                                     
                            Energetic_Disorder_Value.value, 
                                     
                            Above_Gap_Value.value , 
                                     
                            Temperature_Widget.value ) ] )
                           
                        for Energetic_Gap in Sample_E_opts }
        
                #----------------------------------------------------------------------------------------------------------
                # Determine Non-Radiative Loss
                #----------------------------------------------------------------------------------------------------------
            
                if Type == 'Radiative Limit':
                
                    V_oc_Losses = { Energetic_Gap : 0 for Energetic_Gap in Sample_E_opts }
            
                else:
        
                    if Type == 'Empirical - Ullbrich et al.':
            
                        V_oc_Losses = { Energetic_Gap : 
                           
                        Empirical_NR_V_oc_Loss_Ullbrich( Energetic_Gap ) for Energetic_Gap in Sample_E_opts }
                
                    else:
                
                        V_oc_Losses = { Energetic_Gap : 
                               
                        Non_Radiative_Open_Circuit_Voltage_Loss( Energetic_Gap , Type )
                               
                        for Energetic_Gap in Sample_E_opts }
                
                #----------------------------------------------------------------------------------------------------------
                # Determine New Light Power
                #----------------------------------------------------------------------------------------------------------

                if Intensity_Type.value == 'Irradiance (W/m2)':
            
                    Scale_Factor = 1
            
                if Intensity_Type.value == 'Irradiance (mW/cm2)':
            
                    Scale_Factor = 0.1
            
                if Intensity_Type.value == 'Illuminance (lx)':
            
                    Scale_Factor = P_lights[ Spectrum_Type ] / Lux_Values[ Spectrum_Type ][ 'V 2-deg' ]
          
                New_Light_Power = Scale_Factor * Intensity_Input.value
        
                #----------------------------------------------------------------------------------------------------------
                # Determine Figures of Merit
                #----------------------------------------------------------------------------------------------------------
        
                Figures_of_Merit = { Energetic_Gap : 
                                
                    Data_Analyser( EQE_Spectra[ Energetic_Gap ],
                              
                        Spectrum_Type, 
                              
                        New_Light_Power,
                                  
                        Temperature_Widget.value, 
                             
                        V_oc_Losses[ Energetic_Gap ] )
                                
                    for Energetic_Gap in Sample_E_opts }
                    
                Compiled_Figures_of_Merit_Local = Figures_of_Merit_Compiler( Sample_E_opts , Figures_of_Merit )
            
                #----------------------------------------------------------------------------------------------------------
                # Determine Maximum PCE
                #----------------------------------------------------------------------------------------------------------
    
                PCEs = Compiled_Figures_of_Merit_Local[ 'PCE' ]
            
                Optimal_Index = where( PCEs == max( PCEs ) )[ 0 ][ 0 ] 
            
                Optimal_Gap = Sample_E_opts[ Optimal_Index ] 
        
                Best_E_opt = Optimal_Gap
            
                # print( Best_E_opt )

            Energetic_Gaps[ Type ] = Best_E_opt
                        
    #----------------------------------------------------------------------------------------------------------------------
    # Determine Intensities
    #----------------------------------------------------------------------------------------------------------------------  
    
    Intensities = logspace( log10( Min_Irradiance.value ) , log10( Max_Irradiance.value ) , N_Irradiances.value )
    
    global Compiled_Figures_of_Merit_vs_Intensity
    
    Compiled_Figures_of_Merit_vs_Intensity = {}
    
    for Type in Calculation_Types:
        
        Energetic_Gap = Energetic_Gaps[ Type ]
                
        if EQE_Simulation_Type.value == 'Step Function':
            
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Step Function Like
            #--------------------------------------------------------------------------------------------------------------
                
            EQE_Spectrum = array( [
                
                Wavelengths ,
                
                SQ_EQE_Simulator( Energies, 
                                                             
                    Energetic_Gap, 
                                                             
                    Below_Gap_Value.value,
                                                            
                    Above_Gap_Value.value ) ] )
        
        if EQE_Simulation_Type.value == 'Urbach Tail (Inorganic/Perovskite)':
        
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Urbach Tail Discontinuous Transition
            #--------------------------------------------------------------------------------------------------------------
                
            if Use_Thermal_Energy_Checkbox.value == True:
                
                    E_U = Temperature_Widget.value * k / e
                    
            else:
                
                E_U = Urbach_Energy_Value.value
                
            EQE_Spectrum = array( [
                
                Wavelengths ,
                
                E_U_Tail_EQE_Simulator( Energies, 
                                                             
                    Energetic_Gap, 
                                                             
                    E_U,
                                                            
                    Above_Gap_Value.value ) ] )
            
        if EQE_Simulation_Type.value == 'Exciton Absorption (Organic)':
        
            #--------------------------------------------------------------------------------------------------------------
            # Simulate Spectra - Excitonic Absorption
            #--------------------------------------------------------------------------------------------------------------
                
            EQE_Spectrum = array( [
                
                Wavelengths ,
                
                SE_EQE_Simulator( Energies , 
                                     
                    Energetic_Gap, 
                                     
                    Energetic_Disorder_Value.value, 
                                     
                    Above_Gap_Value.value , 
                                     
                    Temperature_Widget.value ) ] )
        
        #--------------------------------------------------------------------------------------------------------------
        # Determine Non-Radiative Loss
        #--------------------------------------------------------------------------------------------------------------
            
        if Type == 'Radiative Limit':
                
            V_oc_Loss = 0
            
        else:
        
            if Type == 'Empirical - Ullbrich et al.':
            
                V_oc_Loss = Empirical_NR_V_oc_Loss_Ullbrich( Energetic_Gap )
                
            else:
                
                V_oc_Loss = Non_Radiative_Open_Circuit_Voltage_Loss( Energetic_Gap , Type )
                               
        #--------------------------------------------------------------------------------------------------------------
        # Determine New Light Power
        #--------------------------------------------------------------------------------------------------------------

        if Intensity_Type.value == 'Irradiance (W/m2)':
            
            Scale_Factor = 1
            
        if Intensity_Type.value == 'Irradiance (mW/cm2)':
            
            Scale_Factor = 0.1
            
        if Intensity_Type.value == 'Illuminance (lx)':
            
            Scale_Factor = P_lights[ Spectrum_Type ] / Lux_Values[ Spectrum_Type ][ 'V 2-deg' ]
          
        Scaled_Intensities = Scale_Factor * Intensities
        
        #--------------------------------------------------------------------------------------------------------------
        # Determine Figures of Merit
        #--------------------------------------------------------------------------------------------------------------
        
        Figures_of_Merit = { Intensities[ i ] : 
                                
            Data_Analyser( EQE_Spectrum,
                              
                Spectrum_Type, 
                              
                Scaled_Intensities[ i ],
                          
                Temperature_Widget.value,
                             
                V_oc_Loss )
                                
            for i in range( len( Intensities ) ) }
                    
        Compiled_Figures_of_Merit_Local = Figures_of_Merit_Compiler_vs_Lux( Intensities , Figures_of_Merit )
            
        Compiled_Figures_of_Merit_vs_Intensity[ Type ] = Compiled_Figures_of_Merit_Local
            
    #--------------------------------------------------------------------------------------------------------------
    # Generate Graphs
    #--------------------------------------------------------------------------------------------------------------
    
    global Output_Graphs_vs_Intensity
    
    Output_Graphs_vs_Intensity = { Curve_Type : Output() for Curve_Type in Curve_Types }
    
    for Curve_Type in Curve_Types:
                
        if Curve_Units[ Curve_Type ] != '':
                    
            Unit = ' (' + Curve_Units[ Curve_Type ] + ')'
                    
        else:
                    
            Unit = ''
                
        with Output_Graphs_vs_Intensity[ Curve_Type ]:
                    
            plt.figure( dpi = 100 )
        
            for Type in Calculation_Types:
                    
                plt.plot( Intensities , Compiled_Figures_of_Merit_vs_Intensity[ Type ][ Curve_Type ] , label = Type )
                                        
            plt.ylabel( Display_Names[ Curve_Type ] + Unit  )
            
            plt.xlabel( Intensity_Type.value )
                    
            if Curve_Type == 'J_0':
                        
                plt.yscale( 'log' )
            
            plt.xscale( 'log' )
            
            plt.legend()
                    
            plt.show()
                    
        Output_Graphs_vs_Intensity[ Curve_Type ]
                    
        Output_Graphs_vs_Intensity[ Curve_Type ].layout.display = 'none'

    #--------------------------------------------------------------------------------------------------------------
    # Output Graphs and Enable Save Data Button
    #--------------------------------------------------------------------------------------------------------------
    
    Save_Intensity_Dep_Data_Button.disabled = False
    
    Copy_Limit_vs_Intensity_Data_to_Clipboard_Button.disabled = False
                    
    Curve_Control = RadioButtons( 
                
        options = Curve_Types,
            
        value = Curve_Types[ 0 ] )
           
    Compiled_Box = VBox( 
        
        [ Curve_Control ] + 
        
        [ Output_Graphs_vs_Intensity[ Curve_Type ] for Curve_Type in Curve_Types ]  )
    
    
    if Determine_Optimal_Gap_Checkbox.value == True:
    
        Loss_Type_Labels = { Loss_Type : Label( value = Loss_Type + ':') for Loss_Type in Calculation_Types }
    
        Best_Energetic_Gap_Labels = { Loss_Type : 
                                 
                                 Label( value = 'E_opt = ' + str( Energetic_Gaps[ Loss_Type ] ) + ' eV' ) 
                                 
                                 for Loss_Type in Calculation_Types }
        
        for Loss_Type in Calculation_Types:
            
            Compiled_Box.children = Compiled_Box.children + ( 
                
                Loss_Type_Labels[ Loss_Type ] , Best_Energetic_Gap_Labels[ Loss_Type ] )

            
    Curve_Control.observe( On_Change_Output_Graph_vs_Intensity_Revealer, names = 'value' )
            
    Curve_Control.value = 'PCE'
                
    Varied_Intensity_Box.children = Varied_Intensity_Box.children[ :-1 ]
            
    Varied_Intensity_Box.children = Varied_Intensity_Box.children + ( Compiled_Box , )

# Similar to the previous function, the intensity-dependent graph that is revealed will be controlled by the following external function: 

def On_Change_Output_Graph_vs_Intensity_Revealer( Change ):
    
    """Change which graph is revealed in the simulation pane.""" 
    
    for Curve_Type in Curve_Types:
        
        Output_Graphs_vs_Intensity[ Curve_Type ].layout.display = 'none'
        
    Output_Graphs_vs_Intensity[ Change[ 'new' ] ].layout.display = None

# 4.7.2. Controlling in the Experimentally-Determined EQE Case

# In this section, the functions that control the simulations in the loaded-in photovoltaic quantum efficiency case are defined. Firstly, two dictionaries are globally defined (such that they can be updated each time the user analyses a spectrum):

Limit_Dep_Figures_of_Merits = {}

Intensity_Dep_Figures_of_Merits = {}

# Following this, a function is defined for conducting simulations with respect to the lower limit of the integral (with respect to photon energy); this is performed by varying the lower limit from data point to data point - no interpolation is performed on the finite data:

def On_Click_Lower_Limit_EQE_Data_Analyser( Button ):
    
    """Analyse the experimental EQE data on click of a button."""
    
    Key = Button.style._view_name
    
    Prior_Graph = Lower_Limit_Graph_Specifiers[ Key ].value

    #-----------------------------------------------------------------------------------------------------------------------
    # Load the Data to Analyse
    #-----------------------------------------------------------------------------------------------------------------------
      
    Wavelengths = EQE_Spectra_to_Investigate[ Key ][ 0 ]
    
    Energies = Energy_Wavelength_Converter( Wavelengths )
    
    EQEs = EQE_Spectra_to_Investigate[ Key ][ 1 ]
    
    EQEs = EQE_Scale_Factors[ Key ].value * EQEs

    #-----------------------------------------------------------------------------------------------------------------------
    # Load Widget Choices
    #-----------------------------------------------------------------------------------------------------------------------
    
    Temperature = Temperature_Widgets[ Key ].value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Loss Type
    #-----------------------------------------------------------------------------------------------------------------------
      
    if Non_Radiative_Loss_Selections[ Key ].value == 'Non-Radiative Loss (V):':
        
        NR_Loss = Voltage_Inputs[ Key ].value
        
    if Non_Radiative_Loss_Selections[ Key ].value == 'Electroluminescent Quantum Efficiency:':
        
        NR_Loss = k * Temperature / e * log( Voltage_Inputs[ Key ].value )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Spectrum Type and Power
    #-----------------------------------------------------------------------------------------------------------------------
    
    Spectrum_Type = Spectrum_Selector.value
    
    if Intensity_Type.value == 'Irradiance (W/m2)':
        
        Scale_Factor = 1 / 10 
        
    if Intensity_Type.value == 'Irradiance (mW/cm2)':
        
        Scale_Factor = 1
                
    if Intensity_Type.value == 'Illuminance (lx)':
        
        V_Type = Luminous_Efficiency_Type_Selector.value
        
        Scale_Factor = 1 / Constants_of_Proportionality[ Spectrum_Type ][ V_Type ] / 10 # Convert to mW/cm2

    Light_Power = Scale_Factor * Intensity_Values[ Key ].value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Photon Fluxes
    #-----------------------------------------------------------------------------------------------------------------------
    
    Photon_Fluxes = Photon_Flux_Interpolator( Spectrum_Type , Wavelengths , Light_Power )
    
    BB_Fluxes = Planck_Photon_Flux_Wavelength( Wavelengths , Temperature )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # If one-sun V_oc has been given, determine figures-of-merit:
    #-----------------------------------------------------------------------------------------------------------------------
    
    if Non_Radiative_Loss_Selections[ Key ].value == 'Experimental Open-Circuit Voltage (V):':
    
        Limit_Dep_Figures_of_Merit = [
        
            One_Sun_V_oc_Data_Analyser( 
                
                array( [ Wavelengths[ :i + 1 ]  , EQEs[ :i + 1 ] ] ), 
            
                Spectrum_Type, 
            
                Light_Power, 
            
                Temperature,
            
                Voltage_Inputs[ Key ].value ) 
    
            for i in range( len( Wavelengths ) ) ]

    #-----------------------------------------------------------------------------------------------------------------------
    # If Non-Radiative Loss or EQE_EL has been given, determine figures-of-merit:
    #-----------------------------------------------------------------------------------------------------------------------
        
    else:
                
        Limit_Dep_Figures_of_Merit = [
        
            Data_Analyser( 
                
                array( [ Wavelengths[ :i + 1 ]  , EQEs[ :i + 1 ] ] ), 
            
                Spectrum_Type, 
            
                Light_Power, 
            
                Temperature,
            
                NR_Loss ) 
    
            for i in range( len( Wavelengths ) ) ]
    
    Current_Curve_Types = list( Limit_Dep_Figures_of_Merit[ 0 ].keys() )     # May Include Delta_V_oc_nr
    
    Limit_Dep_Figures_of_Merit = { Curve_Type :
                                 
        [ Limit_Dep_Figures_of_Merit[ i ][ Curve_Type ] for i in range( len( Wavelengths ) ) ] 
                                 
        for Curve_Type in Current_Curve_Types }
    
    Limit_Dep_Figures_of_Merits[ Key ] = { 'E_lower' : list( Energies ) } | Limit_Dep_Figures_of_Merit
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Graphs
    #-----------------------------------------------------------------------------------------------------------------------
    
    Lower_Limit_Graphs = { Graph_Key : Output() for Graph_Key in Current_Curve_Types }
        
    for Graph_Key in Current_Curve_Types:
        
        with Lower_Limit_Graphs[ Graph_Key ]:
            
            plt.plot( Energies , Limit_Dep_Figures_of_Merit[ Graph_Key ] )
            
            plt.xlabel( 'Lower Limit of Integral, $E_\mathrm{lower}$ (eV)' )
            
            plt.ylabel( Full_Curve_Labels[ Graph_Key ] )
            
            plt.show()
            
    for Graph_Key in Current_Curve_Types:
        
        Lower_Limit_Graphs[ Graph_Key ].layout.display = 'none'
        
    Lower_Limit_FoM_Graphs[ Key ] = Lower_Limit_Graphs
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add Graphs to Box
    #-----------------------------------------------------------------------------------------------------------------------
    
    Graph_Selector = Lower_Limit_Graph_Specifiers[ Key ]
    
    Graph_Selector.observe( Lower_Limit_EQE_Figure_of_Merit_Graph_Selection , names = 'value' )
    
    Graph_Selector.options = Graph_Selector.options[ :1 ]
    
    for Curve_Type in Curve_Types:
        
        Graph_Selector.options = Graph_Selector.options + ( Curve_Type, )
        
    if Non_Radiative_Loss_Selections[ Key ].value == 'Experimental Open-Circuit Voltage (V):':
    
        Graph_Selector.options = Graph_Selector.options[ :3 ] + ( 'Delta_V_oc_nr', ) + Graph_Selector.options[ 3: ]
        
    if Prior_Graph in Current_Curve_Types:
        
        Graph_Selector.value = Prior_Graph
        
    Lower_Limit_Graph_Boxes[ Key ].children = Lower_Limit_Graph_Boxes[ Key ].children[ :1 ] + tuple( [ Lower_Limit_Graphs[ Key ] for Key in Current_Curve_Types ] )
    
    Copy_Data_Buttons[ Key ].disabled = False

# Similarly, at a given lower limit of the integral, the following function will conduct simulations with respect to the intensity of the incident light:

def On_Click_Intensity_EQE_Data_Analyser( Button ):
    
    """Analyse the experimental EQE data on click of a button."""
    
    Key = Button.style._view_name
    
    Prior_Graph = Vs_Intensity_Graph_Specifiers[ Key ].value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Load the Data to Analyse
    #-----------------------------------------------------------------------------------------------------------------------
      
    Wavelengths = EQE_Spectra_to_Investigate[ Key ][ 0 ]
    
    Energies = Energy_Wavelength_Converter( Wavelengths )
    
    EQEs = EQE_Spectra_to_Investigate[ Key ][ 1 ]
    
    EQEs = EQE_Scale_Factors[ Key ].value * EQEs
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Trim the Spectra
    #-----------------------------------------------------------------------------------------------------------------------
          
    Lower_Energy = Vs_Intensity_Lower_Limit_Values[ Key ].value
    
    Good_Indices = [ i for i in range( len( Energies ) ) if Energies[ i ] > Lower_Energy ]
    
    Trimmed_Wavelengths = Wavelengths[ Good_Indices ]
    
    Trimmed_Energies = Energies[ Good_Indices ]
    
    Trimmed_EQEs = EQEs[ Good_Indices ]

    #-----------------------------------------------------------------------------------------------------------------------
    # Load Widget Choices
    #-----------------------------------------------------------------------------------------------------------------------
    
    Temperature = Temperature_Widgets[ Key ].value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Loss Type
    #-----------------------------------------------------------------------------------------------------------------------
      
    if Non_Radiative_Loss_Selections[ Key ].value == 'Non-Radiative Loss (V):':
        
        NR_Loss = Voltage_Inputs[ Key ].value
        
    if Non_Radiative_Loss_Selections[ Key ].value == 'Electroluminescent Quantum Efficiency:':
        
        NR_Loss = k * Temperature / e * log( Voltage_Inputs[ Key ].value )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Spectrum Type and Power
    #-----------------------------------------------------------------------------------------------------------------------
    
    Spectrum_Type = Spectrum_Selector.value
    
    if Intensity_Type.value == 'Irradiance (W/m2)':
        
        Scale_Factor = 1 / 10 
        
    if Intensity_Type.value == 'Irradiance (mW/cm2)':
        
        Scale_Factor = 1
                
    if Intensity_Type.value == 'Illuminance (lx)':
        
        V_Type = Luminous_Efficiency_Type_Selector.value
        
        Scale_Factor = 1 / Constants_of_Proportionality[ Spectrum_Type ][ V_Type ] / 10 # Convert to mW/cm2
      
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Light Powers
    #-----------------------------------------------------------------------------------------------------------------------
    
    Intensities_to_Investigate = logspace( 
    
        log10( Min_Irradiance.value ),
    
        log10( Max_Irradiance.value ), 
        
        N_Irradiances.value )
    
    Correct_Units_Intensities = Scale_Factor * Intensities_to_Investigate
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Determine Photon Fluxes
    #-----------------------------------------------------------------------------------------------------------------------
    
    Photon_Fluxes = [ Photon_Flux_Interpolator( Spectrum_Type , Wavelengths , Light_Power )
                     
        for Light_Power in Correct_Units_Intensities ]
    
    BB_Fluxes = Planck_Photon_Flux_Wavelength( Wavelengths , Temperature )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # If one-sun V_oc has been given, determine figures-of-merit:
    #-----------------------------------------------------------------------------------------------------------------------
    
    if Non_Radiative_Loss_Selections[ Key ].value == 'Experimental Open-Circuit Voltage (V):':
    
        Intensity_Dep_Figures_of_Merit = [
        
            One_Sun_V_oc_Data_Analyser( 
                
                array( [ Trimmed_Wavelengths  , Trimmed_EQEs ] ), 
            
                Spectrum_Type, 
            
                Light_Power, 
            
                Temperature,
            
                Voltage_Inputs[ Key ].value ) 
    
            for Light_Power in Correct_Units_Intensities ]

    #-----------------------------------------------------------------------------------------------------------------------
    # If Non-Radiative Loss or EQE_EL has been given, determine figures-of-merit:
    #-----------------------------------------------------------------------------------------------------------------------
        
    else:
                
        Intensity_Dep_Figures_of_Merit = [
        
            Data_Analyser( 
                
                array( [ Trimmed_Wavelengths  , Trimmed_EQEs ] ), 
            
                Spectrum_Type, 
            
                Light_Power, 
            
                Temperature,
            
                NR_Loss ) 
    
            for Light_Power in Correct_Units_Intensities ]
    
    Current_Curve_Types = list( Intensity_Dep_Figures_of_Merit[ 0 ].keys() )     # May Include Delta_V_oc_nr
    
    Intensity_Dep_Figures_of_Merit = { Curve_Type :
                                 
        [ Intensity_Dep_Figures_of_Merit[ i ][ Curve_Type ] for i in range( len( Correct_Units_Intensities ) ) ] 
                                 
        for Curve_Type in Current_Curve_Types }
    
    Intensity_Dep_Figures_of_Merits[ Key ] = { Intensity_Type.value : list( Intensities_to_Investigate ) } | Intensity_Dep_Figures_of_Merit
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Graphs
    #-----------------------------------------------------------------------------------------------------------------------
    
    Intensity_Graphs = { Graph_Key : Output() for Graph_Key in Current_Curve_Types }
        
    for Graph_Key in Current_Curve_Types:
        
        with Intensity_Graphs[ Graph_Key ]:
            
            plt.plot( Intensities_to_Investigate , Intensity_Dep_Figures_of_Merit[ Graph_Key ] )
            
            plt.xlabel( Intensity_Type.value )
            
            plt.ylabel( Full_Curve_Labels[ Graph_Key ] )
            
            plt.xscale( 'log' )
            
            plt.show()
            
    for Graph_Key in Current_Curve_Types:
        
        Intensity_Graphs[ Graph_Key ].layout.display = 'none'
        
    Intensity_FoM_Graphs[ Key ] = Intensity_Graphs
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add Graphs to Box
    #-----------------------------------------------------------------------------------------------------------------------
    
    Graph_Selector = Vs_Intensity_Graph_Specifiers[ Key ]

    Graph_Selector.observe( Intensity_EQE_Figure_of_Merit_Graph_Selection , names = 'value' )
    
    Graph_Selector.options = Graph_Selector.options[ :1 ]
    
    for Curve_Type in Curve_Types:
        
        Graph_Selector.options = Graph_Selector.options + ( Curve_Type, )
        
    if Non_Radiative_Loss_Selections[ Key ].value == 'Experimental Open-Circuit Voltage (V):':
    
        Graph_Selector.options = Graph_Selector.options[ :3 ] + ( 'Delta_V_oc_nr', ) + Graph_Selector.options[ 3: ]
        
    if Prior_Graph in Current_Curve_Types:
        
        Graph_Selector.value = Prior_Graph
        
    Vs_Intensity_Graph_Boxes[ Key ].children = Vs_Intensity_Graph_Boxes[ Key ].children[ :1 ] + tuple( [ Intensity_Graphs[ Key ] for Key in Current_Curve_Types ] )
    
    Copy_Data_Button_vs_Is[ Key ].disabled = False

# 4.7.3. Spectrum Analysing Tab Generator

# The final control function generates an EQE_PV spectrum-analysing tab each time the user loads in a spectrum.

def EQE_Spectrum_Analysing_Tab_Generator( EQE_Spectrum_Label ):
    
    """Generate the tab for analysing a given EQE spectrum."""
    
    EQE_Spectrum_Analysing_Tab.layout.display = None
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Load the Data to Generate the Plot 
    #-----------------------------------------------------------------------------------------------------------------------
      
    Wavelengths = EQE_Spectra_to_Investigate[ EQE_Spectrum_Label ][ 0 ]
    
    Energies = Energy_Wavelength_Converter( Wavelengths )
    
    EQEs = EQE_Spectra_to_Investigate[ EQE_Spectrum_Label ][ 1 ]

    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Widget for Specifying Whether the Input Voltage is Experimetnal V_oc, or non-rad V_oc loss, or EQE_EL
    #-----------------------------------------------------------------------------------------------------------------------
    
    Non_Radiative_Loss_Selection = RadioButtons( 
        
        options = Non_Rad_Loss_Calc_Types,
        
        value = Non_Rad_Loss_Calc_Types[ 0 ],
    
        style = { '_view_name' : EQE_Spectrum_Label } )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add to the List
    #-----------------------------------------------------------------------------------------------------------------------
    
    Non_Radiative_Loss_Selections[ EQE_Spectrum_Label ] = Non_Radiative_Loss_Selection
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Tell the Widget to Observe the Function for Changing What the User Input is
    #-----------------------------------------------------------------------------------------------------------------------
    
    Non_Radiative_Loss_Selection.observe( Voltage_Input_Changer , names = 'value' )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Widget for Specifying the "Voltage" Input (Can Be EQE_EL depending on above )
    #-----------------------------------------------------------------------------------------------------------------------
    
    Voltage_Input = FloatText(
    
        value = 0.2,
    
        description = 'Non-Radiative Voltage Loss (V):',
        
        style = { 'description_width' : 'initial' } )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add it to the List
    #-----------------------------------------------------------------------------------------------------------------------
    
    Voltage_Inputs[ EQE_Spectrum_Label ] = Voltage_Input
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Widget for Specifying the Temperature and Add it to the List
    #-----------------------------------------------------------------------------------------------------------------------
    
    Temperature_Widget = FloatText( 
    
        value = 293.15 ,
    
        description = 'Temperature (K):',
    
        style = { 'description_width' : 'initial' } )
    
    Temperature_Widgets[ EQE_Spectrum_Label ] = Temperature_Widget
    
    #-----------------------------------------------------------------------------------------------------------------------
    # EQE Scale factor
    #-----------------------------------------------------------------------------------------------------------------------

    EQE_Scale_Factor = FloatText( 
    
        value = 1,
    
        description = 'EQE Scale Factor:',
    
        style = { 'description_width' : 'initial' } )
    
    EQE_Scale_Factors[ EQE_Spectrum_Label ] = EQE_Scale_Factor
    
    #-----------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------
    # Versus Lower Limit Widgets
    #-----------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Widget for Selecting Which Graph is Shown (Figure-of-Merit Graphs Will be Added Later)
    #-----------------------------------------------------------------------------------------------------------------------

    Lower_Limit_Graph_Specifier = RadioButtons( 
    
        options = [ 'EQE' ],
    
        value = 'EQE',
    
        style = { '_view_name' : EQE_Spectrum_Label } )
    
    Lower_Limit_Graph_Specifiers[ EQE_Spectrum_Label ] = Lower_Limit_Graph_Specifier
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Plot the EQE Spectrum in a Graph and Add Graph to List
    #-----------------------------------------------------------------------------------------------------------------------

    Lower_Limit_EQE_Graph = Output()
    
    with Lower_Limit_EQE_Graph:
        
        plt.plot( Energies , EQEs )
        
        plt.xlabel( 'Photon Energy, E (eV)' )
        
        plt.ylabel( 'EQE' )
        
        plt.show()
    
    Lower_Limit_EQE_Graphs[ EQE_Spectrum_Label ] = Lower_Limit_EQE_Graph 
    
    Vs_Intensity_EQE_Graph = Output()
    
    with Vs_Intensity_EQE_Graph:
        
        plt.plot( Energies , EQEs )
        
        plt.xlabel( 'Photon Energy, E (eV)' )
        
        plt.ylabel( 'EQE' )
        
        plt.show()
    
    Vs_Intensity_EQE_Graphs[ EQE_Spectrum_Label ] = Vs_Intensity_EQE_Graph 
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add Graph to Box, Add Box to List of Graph Boxes
    #-----------------------------------------------------------------------------------------------------------------------
    
    Lower_Limit_Graph_Box = VBox( [ Lower_Limit_EQE_Graph ] )
    
    Lower_Limit_Graph_Boxes[ EQE_Spectrum_Label ] = Lower_Limit_Graph_Box
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Two Panels - First One With Overall Input Parameters , the Other with Different Simulation Types
    #-----------------------------------------------------------------------------------------------------------------------
    
    Left_Hand_Box = Tab()
    
    Left_Hand_Box.children = [
        
        VBox( [ 
            
            Spectrum_Selector_Box,
        
            Temperature_Widget,
        
            EQE_Scale_Factor,
        
            Non_Radiative_Loss_Selection,
        
            Voltage_Input
        
        ] ) ]
    
    Left_Hand_Box.set_title( 0 , 'Input Parameters' )
        
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Two Panels - Second One with Different Simulation Types (vs lower limit and vs intensity)
    #-----------------------------------------------------------------------------------------------------------------------
  
    Intensity_Value = FloatText( 
        
        value = 100,
    
        description = Intensity_Type.value,
    
        style = { 'description_width' : 'initial' } )
    
    Intensity_Values[ EQE_Spectrum_Label ] = Intensity_Value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Button for Analysing the Data and Add it to the List
    #-----------------------------------------------------------------------------------------------------------------------

    Analyse_Lower_Limit_Data_Button = Button( 
        
        description = 'Analyse Data' , 
    
        style = { '_view_name' : EQE_Spectrum_Label } )
    
    Analyse_Lower_Limit_Data_Button.on_click( On_Click_Lower_Limit_EQE_Data_Analyser )

    Analyse_Lower_Limit_Data_Buttons[ EQE_Spectrum_Label ] = Analyse_Lower_Limit_Data_Button
    
    Copy_Data_Button = Button( 
        
        description = 'Copy to Clipboard',
        
        style = { '_view_name' : EQE_Spectrum_Label },
    
        disabled = True )
    
    Copy_Data_Button.on_click( On_Click_EQE_Versus_E_Lower_to_Clipboard_Copier )
    
    Copy_Data_Buttons[ EQE_Spectrum_Label ] = Copy_Data_Button
        
    Right_Hand_Box_vs_E_lower = VBox( [ 
            
        Intensity_Type,
            
        Intensity_Value,
            
        HBox( [ Analyse_Lower_Limit_Data_Button, Copy_Data_Button ] ),
            
        Lower_Limit_Graph_Specifier,
            
        Lower_Limit_Graph_Box ] )
    
    #-----------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------
    # Versus Intensity Widgets
    #-----------------------------------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------------------------------
    
    Vs_Intensity_Lower_Limit_Value = FloatText(
    
        value = 1,
    
        description = 'Lower Limit of Integral, $E_\mathrm{lower}$ (eV):',
        
        style = { 'description_width' : 'initial' } )
    
    Vs_Intensity_Lower_Limit_Values[ EQE_Spectrum_Label ] = Vs_Intensity_Lower_Limit_Value
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Widget for Selecting Which Graph is Shown (Figure-of-Merit Graphs Will be Added Later)
    #-----------------------------------------------------------------------------------------------------------------------

    Vs_Intensity_Graph_Specifier = RadioButtons( 
    
        options = [ 'EQE' ],
    
        value = 'EQE',
    
        style = { '_view_name' : EQE_Spectrum_Label } )
    
    Vs_Intensity_Graph_Specifiers[ EQE_Spectrum_Label ] = Vs_Intensity_Graph_Specifier
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Generate Graph Box
    #-----------------------------------------------------------------------------------------------------------------------
    
    Vs_Intensity_Graph_Box = VBox( [ Vs_Intensity_EQE_Graph ] )
    
    Vs_Intensity_Graph_Boxes[ EQE_Spectrum_Label ] = Vs_Intensity_Graph_Box
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create the Button for Analysing the Data and Add it to the List
    #-----------------------------------------------------------------------------------------------------------------------

    Analyse_Intensity_Data_Button = Button( 
        
        description = 'Analyse Data' , 
    
        style = { '_view_name' : EQE_Spectrum_Label } )
    
    Analyse_Intensity_Data_Button.on_click( On_Click_Intensity_EQE_Data_Analyser )

    Analyse_Intensity_Data_Buttons[ EQE_Spectrum_Label ] = Analyse_Intensity_Data_Button
    
    Copy_Data_Button_vs_I = Button(
        
        description = 'Copy to Clipboard',
                                   
        style = { '_view_name' : EQE_Spectrum_Label },
    
        disabled = True )
    
    Copy_Data_Button_vs_I.on_click( On_Click_EQE_Versus_Intensity_to_Clipboard_Copier )
    
    Copy_Data_Button_vs_Is[ EQE_Spectrum_Label ] = Copy_Data_Button_vs_I
    
    Right_Hand_Box_vs_Intensity = VBox( [ 
            
        Vs_Intensity_Lower_Limit_Value,
            
        Intensity_Type,
            
        Min_Irradiance, 
            
        Max_Irradiance, 
            
        N_Irradiances,
            
        HBox( [ Analyse_Intensity_Data_Button , Copy_Data_Button_vs_I ] ),
            
        Vs_Intensity_Graph_Specifier,
            
        Vs_Intensity_Graph_Box ] ) 
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Create Tab
    #-----------------------------------------------------------------------------------------------------------------------

    Right_Hand_Box = Tab()
    
    Right_Hand_Box.children = (
        
        Right_Hand_Box_vs_E_lower,
        
        Right_Hand_Box_vs_Intensity )
    
    Right_Hand_Box.set_title( 0 , 'Versus Lower Limit' )

    Right_Hand_Box.set_title( 1 , 'Versus Intensity' )
    
    Output_Box = HBox( [ Left_Hand_Box , Right_Hand_Box ] )
    
    #-----------------------------------------------------------------------------------------------------------------------
    # Add Box to Analysis Tab
    #-----------------------------------------------------------------------------------------------------------------------
    
    EQE_Spectrum_Analysing_Tab.children = EQE_Spectrum_Analysing_Tab.children + ( Output_Box , )
    
    N_Contents = len( EQE_Spectrum_Analysing_Tab.children )
    
    EQE_Spectrum_Analysing_Tab.set_title( N_Contents - 1 , EQE_Spectrum_Label )

# 5. The User Interface

# The user interface is essentially broken into two main parts - the simulated photovoltaic external quantum efficiency spectrum-analysing part and the experimental one; the former is then broken into three further parts. The first of these are overall inputs like which photon energies should be used to simulate the EQE spectra, the second are inputs for the optical gap-dependent simulations, and the third are the inputs for the irradiance/illuminance-dependent simulations.

# 5.1. Simulated EQE Interface

# The power conversion efficiency-simulating interface is broken into three parts - the first part includes the inputs for the simulation, the second simulates power conversion efficiency versus optical gap at a particular lux value, whereas the other simulates power conversion efficiency versus lux at a particular optical gap. Both do their respective jobs using one of the particular models outlined above.

# 5.1.1. Overall Inputs

# The overall inputs for the power conversion efficiency-simulating tool are (i) the incident spectrum type, (ii) the intensity unit, (iii) the temperature,(iv) the model for the sub-gap photovoltaic quantum efficiency and (v) the photon energies to simulate it at, and (vi) the non-radiative loss type. The widget for selecting the spectrum was defined as "Spectrum_Selector" in the previous section, as was the intensity unit-selecting widget (as "Intensity_Type") and the temperature-specifying widget. The first new collection of widgets defined in this section are those needed to simulate the photovoltaic external quantum efficiency spectrum and photon energies. The minimum photon energy widget is specified as:

Min_Sim_Energy = FloatText( 
    
    value = 0.001 ,
    
    min = 0.01,
                          
    description = 'Minimum Photon Energy (eV):',

    style = { 'description_width' : 'initial' } )

# The maximum photon energy widget is created uing:

Max_Sim_Energy = FloatText( 
    
    min = 0.01,
    
    value = 5 ,
                          
    description = 'Maximum Photon Energy (eV):',

    style = { 'description_width' : 'initial' } )

# Following this, the widget for customising the number of points to simulate at is created using:

N_Sim_Energies = IntText( 
    
    value = 1001 ,
    
    min = 5,
                          
    description = 'Number of Photon Energies:',

    style = { 'description_width' : 'initial' } )

# These three photon energy-specifying widgets are then compiled into one box:

Simulation_Energies_Box = VBox( [
    
    Min_Sim_Energy,
    
    Max_Sim_Energy,
    
    N_Sim_Energies ] )

# Following this, a function is defined to simulate the photon energies using the above widget values.

def Simulation_Wavelengths_Maker( Min_Sim_Energy , Max_Sim_Energy, N_points ):
    
    """Create linearly-spaced wavelengths in nanometres."""
    
    Energies = linspace( Min_Sim_Energy , Max_Sim_Energy , N_points )
    
    return Energy_Wavelength_Converter( Energies )[ ::-1 ]

# A radio button widget is now defined for specifying whether the absorption is being modelled as a step function (Shockley-Queisser mode), a sub-gap "Urbach tail" model, or the exciton absorption model:

EQE_Simulation_Type = RadioButtons( 
    
    options = [ 'Step Function', 'Urbach Tail (Inorganic/Perovskite)', 'Exciton Absorption (Organic)' ],

    description = 'EQE Model:',

    style = { 'description_width' : 'initial' } ) 

# In the most rudimental of the models for the photovoltaic external quantum efficiency, a step function is used. The above-gap value for this step function is specified as:

Above_Gap_Value = FloatText( 
    
    min = 0,
    
    max = 1,

    value = 0.8,

    description = 'Above-Gap EQE:',

    style = { 'description_width' : 'initial' } )

# Whereas the below-gap value is specified using:

Below_Gap_Value = FloatText( 

    value = 0,
    
    min = 0,
    
    max = 1,
    
    description = 'Below-Gap EQE:',

    style = { 'description_width' : 'initial' } )

# On the other, if a sub-gap "Urbach tail" model is used, the corresponding Urbach energy of the tail can be specified using:

Urbach_Energy_Value = FloatText( 

    value = 0.05,
    
    min = 0,
    
    max = 1,
    
    description = 'Urbach Energy (eV):',

    style = { 'description_width' : 'initial' } )

# As written in the manuscript, previous works by the authors suggest that the Urbach energy is equal to the thermal energy in organic semiconductors. For this reason, the following checkbox is defined to set the Urbach energy equal to the thermal energy <br/><br/> 

Use_Thermal_Energy_Checkbox = Checkbox(

    value = False,

    description = 'Use Thermal Energy' )

# Initially, the script assumes that the User wants to use a step function-type EQE. These widgets for controlling the the Urbach tail parameters are therefore initially hidden using:

Use_Thermal_Energy_Checkbox.layout.display = 'none'

Urbach_Energy_Value.layout.display = 'none'

# The following function is defined such that, if the the User chooses to set the Urbach enegy equal to the thermal energy, the Urbach energy input widget is disabled:

def On_Change_Use_Thermal_Energy( Change ):
    
    """Disable the option to customise the Urbach energy if the thermal energy is to be used."""
    
    if Change[ 'new' ] == True:
        
        Urbach_Energy_Value.disabled = True
        
    if Change[ 'new' ] == False:
        
        Urbach_Energy_Value.disabled = False

# The checkbox widget is then told to obey this function using:

Use_Thermal_Energy_Checkbox.observe( On_Change_Use_Thermal_Energy , names = 'value' )

# With the widgets for customising the sub-gap Urbach tail complete, one last widgets is defined for customising the energetic disorder associated with the sub-gap EQE in the excitonic absorption model:
 
Energetic_Disorder_Value = FloatText( 

    value = 0.05,
    
    min = 0,
    
    description = 'Excitonic Static Disorder (eV):',

    style = { 'description_width' : 'initial' } )

Energetic_Disorder_Value.layout.display = 'none'

# With all the widgets for customising the EQE now specified, they are compiled into a single box:

EQE_Simulating_Box = VBox( [
    
    EQE_Simulation_Type,
    
    Above_Gap_Value,
    
    Below_Gap_Value,

    Urbach_Energy_Value , 
    
    Use_Thermal_Energy_Checkbox,

    Energetic_Disorder_Value ] )

# The following function is used to change which input boxes are revealed, depending on which EQE model the user has selected: <br/><br/> 

def On_Change_EQE_Input_Box_Hider( Change ):
    
    """Hide/reveal the input widgets according to the type of sub-gap absorption."""
    
    New_Option = Change[ 'new' ]
    
    if New_Option == 'Step Function':
                
        Below_Gap_Value.layout.display = None
        
        Urbach_Energy_Value.layout.display = 'none'
        
        Use_Thermal_Energy_Checkbox.layout.display = 'none'

        Energetic_Disorder_Value.layout.display = 'none'
        
    if New_Option == 'Urbach Tail (Inorganic/Perovskite)':
        
        Below_Gap_Value.layout.display = 'none'
        
        Urbach_Energy_Value.layout.display = None
        
        Use_Thermal_Energy_Checkbox.layout.display = None
        
        Energetic_Disorder_Value.layout.display = 'none'
        
    if New_Option == 'Exciton Absorption (Organic)':
        
        Below_Gap_Value.layout.display = 'none'
        
        Urbach_Energy_Value.layout.display = 'none'
        
        Use_Thermal_Energy_Checkbox.layout.display = 'none'
        
        Energetic_Disorder_Value.layout.display = None

EQE_Simulation_Type.observe( On_Change_EQE_Input_Box_Hider , names = 'value' )

# With the widgets for simulating the photovoltaic external quantum efficiency defined, a pseudo-widget for controlling the non-radiative open-circuit voltage calculation type is defined as a list of checkbox widgets. This allows multiple types of open-circuit voltage losses to be considered in a single plot:
 
PCE_Calculation_Types = [ 'Radiative Limit' ] +  Non_Radiative_Loss_Types + [ 'Empirical - Ullbrich et al.']

PCE_Calculation_Select_Widgets = { Type : 
                                 
    Checkbox( value = False,
            
            description = Type,
            
            style = { 'description_width' : 'initial' ,
                    
                    '_view_name' : Type } )
                                 
    for Type in PCE_Calculation_Types }

PCE_Calculation_Selection_Box = VBox( [
    
    Label( value = 'Calculation Type:' ) ] +
    
    [ PCE_Calculation_Select_Widgets[ Type ] for Type in PCE_Calculation_Types ] )

# Initially, only the radiative limit is assumed to be wanted, so the corresponding checkbox's value is set as 'True':

PCE_Calculation_Select_Widgets[ 'Radiative Limit' ].value = True 

# These overall inputs are now compiled into a single 'tab' widget using:

Overall_Inputs = Tab()

Overall_Inputs.children = tuple( [
    
    VBox( [ 
    
        Spectrum_Selector_Box,
    
        Temperature_Widget, ] ),
    
    Simulation_Energies_Box,
    
    EQE_Simulating_Box,
    
    PCE_Calculation_Selection_Box ] )

# The tab titles are specified using:
 
Overall_Inputs_Titles = [ 'Spectrum' , 'Energies' , 'EQEs' , 'NR Loss' ]

for i in range( len( Overall_Inputs_Titles ) ):
    
    Overall_Inputs.set_title( i , Overall_Inputs_Titles[ i ] )

# 5.1.2. Power Conversion Efficiency-Simulating Interface for Varied Optical Gaps

# With the overall inputs now specified, widgets can be created for simulating as a function of the optical gap. Firstly, a widget is created for specifying the total irradiance/illuminance value of the simulation:

Intensity_Input = FloatText(

    value = 100,

    description = 'Intensity (mW/cm2):',

    style = { 'description_width' : 'initial' } )

# Next, three widgets are defined for specifying the minimum optimal gap, the maximum optical gap, and the number of optical gaps to simulate photovoltaic external quantum efficiency (and PCE values) for:

Min_Energetic_Gap = FloatText(
    
    value = 0.7,

    description = 'Min. Energetic Gap (eV)',

    style = { 'description_width' : 'initial' } )

Max_Energetic_Gap = FloatText(
    
    value = 3,

    description = 'Max. Energetic Gap (eV)',

    style = { 'description_width' : 'initial' } )

N_Energetic_Gaps = IntText(
    
    value = 101,

    description = 'Number of Energetic Gaps',

    style = { 'description_width' : 'initial' } )

# These values are compiled into a box to make calling all of them at once a bit easier:

Energetic_Gaps = VBox( [ 
    
    Min_Energetic_Gap,
    
    Max_Energetic_Gap,

    N_Energetic_Gaps] )

# All the widgets needed to simulate figures-of-merit as a function of bandgap are now compiled into a single box (with a placeholder for graphs and three buttons - one for calculating, one for saving data, and the last for copying data to clipboard):

Calculate_Limit_Button = Button( description = 'Compute Limits' )

Calculate_Limit_Button.on_click( On_Click_Limit_vs_E_opt_Computer )

Save_E_opt_Dep_Data_Button = Button( description = 'Save Data' , disabled = True )

Save_E_opt_Dep_Data_Button.on_click( On_Click_Figure_of_Merit_vs_E_opt_Saver )

Copy_Limit_vs_E_opt_Data_to_Clipboard_Button = Button( description = 'Copy to Clipboard' , disabled = True )

Copy_Limit_vs_E_opt_Data_to_Clipboard_Button.on_click( On_Click_E_opt_Dep_Data_to_Clipboard_Copier )

# All these widgets are then compiled into one box:

Varied_E_opt_Box = VBox( [ 
    
    Intensity_Type , 
    
    Intensity_Input, 
    
    Placeholder_Label,
    
    Energetic_Gaps , 

    HBox( [ Calculate_Limit_Button , Copy_Limit_vs_E_opt_Data_to_Clipboard_Button , Save_E_opt_Dep_Data_Button ] ),

    Placeholder_Label ] )

# 5.1.3. Power Conversion Efficiency-Simulating Interface for Varied Irradiances

# Similar to above, widgets are defined for minimum and maximum values, and number of points, but for intensity rather than optical gap:

Min_Irradiance = FloatText( 
    
    value = 1 ,
    
    min = 1e-40,
                          
    description = 'Minimum Irradiance (W/m2):',

    style = { 'description_width' : 'initial' } )

Max_Irradiance = FloatText( 
    
    value = 10000 ,
                          
    description = 'Maximum Irradiance (W/m2):',

    style = { 'description_width' : 'initial' } )

N_Irradiances = IntText( 
    
    value = 201 ,
    
    min = 5,
                          
    description = 'Number of Irradiances:',

    style = { 'description_width' : 'initial' } )

# These widgets are also compiled into an irradiance controlling box:

Irradiances_Box = VBox( [
    
    Min_Irradiance,
    
    Max_Irradiance,
    
    N_Irradiances ] )

# In addition, an optical gap specifying widget is created:

E_opt_Value = FloatText( 
    
    value = 1.5 ,
    
    min = 0,
                          
    description = 'Optical Gap (eV):',

    style = { 'description_width' : 'initial' } )

# One final checkbox widget is defined, which instructs the code to find the optimal optical gap, and plot PCE versus lux for this optimal gap.

Determine_Optimal_Gap_Checkbox = Checkbox( 
    
    value = False,

    description = 'Find and Use Best Gap',

    style = { 'description_width' : 'initial' } )

# If the User chooses for the code to find the best optical gap (that gives the highest PCE), the option to enter a cutomised value is disabled:

def On_Change_Find_Best_Gap( Change ):
    
    """Disable the option to customise the Urbach energy if the thermal energy is to be used."""
    
    if Change[ 'new' ] == True:
        
        E_opt_Value.description = 'Decimal Precision:'
        
        E_opt_Value.value = 5
        
    if Change[ 'new' ] == False:
        
        E_opt_Value.description = 'Optical Gap (eV):'
        
        E_opt_Value.value = 1.5

# The above function is implemented using:

Determine_Optimal_Gap_Checkbox.observe( On_Change_Find_Best_Gap , names = 'value' ) 

# Three more buttons are now defined to finish this section, one for carrying out the calculations, one for saving data, and the other for copying the data quickly to the clipboard:

Calculate_Limit_vs_Intensity_Button = Button( description = 'Compute Limits' )

Calculate_Limit_vs_Intensity_Button.on_click( On_Click_Limit_vs_Intensity_Computer )

Save_Intensity_Dep_Data_Button = Button( description = 'Save Data' , disabled = True )

Save_Intensity_Dep_Data_Button.on_click( On_Click_Figure_of_Merit_vs_Intensity_Saver )

Copy_Limit_vs_Intensity_Data_to_Clipboard_Button = Button( description = 'Copy to Clipboard' , disabled = True )

Copy_Limit_vs_Intensity_Data_to_Clipboard_Button.on_click( On_Click_Intensity_Dep_Data_to_Clipboard_Copier )

# All these widgets are then compiled into a single box:

Varied_Intensity_Box = VBox( [
    
    E_opt_Value,
    
    Determine_Optimal_Gap_Checkbox,
    
    Placeholder_Label,
    
    Intensity_Type,
    
    Irradiances_Box,
    
    HBox( [ 
        
        Calculate_Limit_vs_Intensity_Button , 
        
        Copy_Limit_vs_Intensity_Data_to_Clipboard_Button ,
        
        Save_Intensity_Dep_Data_Button ] ),
    
    Placeholder_Label ] )

# 5.2. EQE-Spectrum Analysing Interface

# In this section, the interface for analysing imported EQE_PV spectra is woven together. This begins by defining the ways to quantify non-radiative open-circuit voltage losses:

Non_Rad_Loss_Calc_Types = [ 'Non-Radiative Loss (V):' , 
                           
                           'Experimental Open-Circuit Voltage (V):',
                          
                           'Electroluminescent Quantum Efficiency:' ]



# The above options will be unique to each imported spectrum, the input value is altered through the following function:

def Voltage_Input_Changer( Change ):
    
    """Change the name of the input voltage parameter."""
    
    Key = Change[ 'owner' ].style._view_name
    
    V_in  = Voltage_Inputs[ Key ]
    
    if Change[ 'new' ] == 'Non-Radiative Loss (V):':
        
        V_in.value  = 0.2
        
        V_in.description = 'Non-Radiative Voltage Loss (V):'
        
    if Change[ 'new' ] == 'Experimental Open-Circuit Voltage (V):':
        
        V_in.value = 1
        
        V_in.description = 'Experimental Open-Circuit Voltage (V):'
        
    if Change[ 'new' ] == 'Electroluminescent Quantum Efficiency:':
        
        V_in.value = 1E-3
        
        V_in.description = 'Electroluminescent Quantum Efficiency:'

# Two functions are now defined for updating the graphs displayed in the spectrum-analysing interface. The first of these changes the graphs in the case of a varied lower limit at fixed intensity: 

def Lower_Limit_EQE_Figure_of_Merit_Graph_Selection( Change ):
    
    """Change which graph is revealed in the EQE-analyser toolbox."""
    
    Key = Change[ 'owner' ].style._view_name
    
    Current_Curve_Types = list( Lower_Limit_Graph_Specifiers[ Key ].options )
    
    Current_Curve_Types.remove( 'EQE' )
    
    if Change[ 'new' ] == 'EQE':
        
        Lower_Limit_EQE_Graphs[ Key ].layout.display = None
        
        for Curve_Type in Current_Curve_Types:
            
            Lower_Limit_FoM_Graphs[ Key ][ Curve_Type ].layout.display == 'none'
            
    else:
            
        Lower_Limit_EQE_Graphs[ Key ].layout.display = 'none'
        
        New_Type = Change[ 'new' ]

        for Curve_Type in Current_Curve_Types:
            
            if Curve_Type != New_Type:
                
                Lower_Limit_FoM_Graphs[ Key ][ Curve_Type ].layout.display = 'none'
                
            else:
                
                Lower_Limit_FoM_Graphs[ Key ][ Curve_Type ].layout.display = None

# Whereas the other changes the graphs in the case of a varied intensity (at fixed optical gap):

def Intensity_EQE_Figure_of_Merit_Graph_Selection( Change ):
    
    """Change which graph is revealed in the EQE-analyser toolbox."""
    
    Key = Change[ 'owner' ].style._view_name
    
    Current_Curve_Types = list( Vs_Intensity_Graph_Specifiers[ Key ].options )
    
    Current_Curve_Types.remove( 'EQE' )
    
    if Change[ 'new' ] == 'EQE':
        
        Vs_Intensity_EQE_Graphs[ Key ].layout.display = None
        
        for Curve_Type in Current_Curve_Types:
            
            Intensity_FoM_Graphs[ Key ][ Curve_Type ].layout.display == 'none'
            
    else:
            
        Vs_Intensity_EQE_Graphs[ Key ].layout.display = 'none'
        
        New_Type = Change[ 'new' ]

        for Curve_Type in Current_Curve_Types:
            
            if Curve_Type != New_Type:
                
                Intensity_FoM_Graphs[ Key ][ Curve_Type ].layout.display = 'none'
                
            else:
                
                Intensity_FoM_Graphs[ Key ][ Curve_Type ].layout.display = None

# Finally, a bunch of dictionaries are defined for storing widgets, the values held by widgets, the graphs, and the compiled boxes. These dictionaries are updated with each additional spectrum that is imported into the interface:

Non_Radiative_Loss_Selections = {}
    
Voltage_Inputs = {}

Temperature_Widgets = {}

Intensity_Values = {}

Analyse_Lower_Limit_Data_Buttons = {}

Lower_Limit_Graph_Specifiers = {}

EQE_Scale_Factors = {}

Lower_Limit_EQE_Graphs = {}

Lower_Limit_FoM_Graphs = {}

Lower_Limit_Graph_Boxes = {} #.append FoM graphs + Hidden

Copy_Data_Buttons = {}

Vs_Intensity_EQE_Graphs = {}

Vs_Intensity_Lower_Limit_Values = {}

Vs_Intensity_Graph_Specifiers = {}

Vs_Intensity_Graph_Boxes = {}

Analyse_Intensity_Data_Buttons = {}

Intensity_FoM_Graphs = {}

Copy_Data_Button_vs_Is = {}

# Finally, an empty experimental EQE_PV spectrum-analysing tab is made (it will be populated each time the user adds a spectrum), then hidden:

EQE_Spectrum_Analysing_Tab = Tab()
EQE_Spectrum_Analysing_Tab.layout.display = 'none'

# 5.3. Compiling Components of the User Interface

# The simulating box is compiled as:

PCE_Calculating_Box = Tab()

PCE_Calculating_Box.children = [
    
    Varied_E_opt_Box,
    
    Varied_Intensity_Box ]

PCE_Calculating_Box.set_title( 0 , 'Varied Optical Gap' )

PCE_Calculating_Box.set_title( 1 , 'Varied Intensity' )

Simulating_Interface = HBox( [ Overall_Inputs , PCE_Calculating_Box ] )

# Update Progress

Progress.value = 100

# Combine this box with the simulating box to make the full interface:

Full_Interface = Tab()

Full_Interface.children = [
    
    Simulating_Interface,
    
    VBox( [ 
        
        Compiled_Add_Spectrum_Box,
    
        EQE_Spectrum_Analysing_Tab ] ) ]

Full_Interface.set_title( 0 , 'Simulation Interface' )

Full_Interface.set_title( 1 , 'Analysis Interface' )

# 5.4. The Interface

Progress.layout.display = 'none'

display( Full_Interface )
