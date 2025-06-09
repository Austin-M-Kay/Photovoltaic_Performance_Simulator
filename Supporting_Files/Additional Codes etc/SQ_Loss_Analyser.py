# -*- coding: utf-8 -*-
"""
Created on Tue Apr 22 15:56:13 2025

@author: austi
"""

def SQ_Loss_Analyser( E_g, Light_Power, Incident_Photon_Fluxes, Black_Body_Fluxes,
    
                     Energies, Wavelengths, kT,
    
                     EQE_EL = 1,
                     
                     R_sh = inf ):
    
    """For a given bandgap (E_g) determine the PCE and the various losses 
    
    in the SQ model."""
    
    EQEs = array( [ 1 if E > E_g else 0 for E in Energies ] )

    J_ph = e * simps( EQEs * Incident_Photon_Fluxes, Wavelengths ) # mA/cm2

    J_0 = e * simps( EQEs * Black_Body_Fluxes, Wavelengths ) / EQE_EL # mA/cm2

    V_oc = kT * log( J_ph / J_0 + 1 )
    
    #--------------------------------------------------------------------------
    # Evaluate PCE
    #--------------------------------------------------------------------------    
    
    def Function_to_minimize( V ):
        
        return V * J_i_Calculator( V, EQEs, Incident_Photon_Fluxes,
            
                                  Black_Body_Fluxes, Energies, Wavelengths,
            
                                  EQE_EL, kT, R_sh,
            
                                  Coupling = False,  # Only considering single junction
                                  
                                  EQEs_em = 1,
            
                                  V_em = 0  )
    
    x0 = 0.8 * V_oc
    
    V_mpp = minimize( Function_to_minimize, x0, bounds = ( ( 0, V_oc ), ) )[ 'x' ][ 0 ]
    
    J_mpp = J_i_Calculator( V_mpp, EQEs, Incident_Photon_Fluxes,
        
                              Black_Body_Fluxes, Energies, Wavelengths,
        
                              EQE_EL, kT, R_sh,
        
                              Coupling = False,  # Only considering single junction
                              
                              EQEs_em = 1,
        
                              V_em = 0  )
    
    P_mpp = abs( J_mpp * V_mpp ) # mW/cm2
    
    FF = P_mpp / J_ph / V_oc

    PCE = P_mpp / Light_Power
    
    #--------------------------------------------------------------------------
    # Evaluate Not Absorbed Light
    #--------------------------------------------------------------------------    
    
    P_unabs = e * simps( Energies * Incident_Photon_Fluxes * ( 1 - EQEs ), Wavelengths )
    
    Frac_unabs = P_unabs / Light_Power
    
    #--------------------------------------------------------------------------
    # Evaluate Thermalised Light
    #--------------------------------------------------------------------------    
        
    Delta_Es = Energies - E_g 
        
    Delta_Es = array( [ DE if DE > 0 else 0 for DE in Delta_Es ] )
    
    P_therm = e * simps( Delta_Es * EQEs * Incident_Photon_Fluxes, Wavelengths )
    
    Frac_therm = P_therm / Light_Power
    
    #--------------------------------------------------------------------------
    # Evaluate Other Losses
    #--------------------------------------------------------------------------    
    
    P_other = Light_Power - P_mpp - P_unabs - P_therm
    
    Frac_other = P_other / Light_Power
    
    return { 'J_sc' : J_ph,
              
               'V_oc' : V_oc,
               
               'V_mpp' : V_mpp,
               
               'J_mpp' : J_mpp,
               
               'P_mpp' : P_mpp,
               
               'FF' : FF,
               
               'PCE' : PCE,
               
               'P_unabs' : P_unabs,
               
               'Frac_unabs': Frac_unabs,
               
               'P_therm' : P_therm,
               
               'Frac_therm' : Frac_therm,
               
               'P_other' : P_other,
               
               'Frac_other' : Frac_other }