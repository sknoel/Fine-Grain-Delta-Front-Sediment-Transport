# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 13:08:16 2021

@author: sknoel
"""

#FUNCTIONS
#General Imports
import math as m

#%%

#REACH SLOPE CALCULATION OF BED SHEAR STRESS
#Input: Depth (m); Reach Slope (unitless)
#Output: Pa
def Tau_DepthSlope(depth, slope):
    
    #Function Constants
    g = 9.81            #gravitationa acceleration (m/s2)
    rho_f = 1.00        #fluid velocity (g/cm3)
    
    BedShearStress = rho_f * g * depth * slope
    
    return BedShearStress

#%%

#DEPTH AVERAGED LAW OF THE WALL FLUID SHEAR STRESS
#Input: Velocity (m/s); Height (m); Saltation Height (m)
#Output: Pa
def Tau_LotW(velocity, height, SaltationHeight):
    
    #Function Constants
    k = 0.4         #von Karmen Constant
    rho_f = 1000    #density of water (kg/m3)
    
    FluidShearVelocity = ( velocity * k * (height - SaltationHeight) ) / ( height * m.log ( height / SaltationHeight ) - height + SaltationHeight)
    FluidShearStress = ( FluidShearVelocity ** 2 ) * rho_f
    
    return FluidShearStress

#%%

#FISCHENICH (2001) CRITICAL SHEAR STRESS
#Input: Diameter (microns)
#Output: Pa
def Tau_CR(diameter):
    
    #Function Constants
    g = 9.81                #gravitational acceleration (m/s2)
    nu = 10**(-6)           #kinematic viscosity of water
    rho_f = 1000            #fluid density (kg/m3)
    rho_s = 2650            #solid density (kg/m3)    
    
    SG = rho_s / rho_f                                              #grain specific gravity (unitless)
    diameter = diameter * 10**(-6)                                  #diameter converstion (m)
    d_star = diameter * ( (SG - 1) * g / nu**2 ) ** (1/3)           

    if diameter > 0.0005:
        phi_r = 31
    
    elif diameter > 0.001:
        phi_r = 32
    
    elif diameter > 0.002:
        print('Gravel Size Sediment, Enter Sand or Smaller')
        
    else:
        phi_r = 30
    
      
    if 2*10**(-6) < diameter <= 4*10**(-6):       #Clays
        Tau_CR = 0.5 * g * (rho_s - rho_f) * diameter * m.tan(phi_r * m.pi/180)
        
    elif 4*10**(-6) < diameter < 0.002:       #Silts & Sands
        Tau_CR = 0.25 * d_star**(-0.6) * g * (rho_s - rho_f) * diameter * m.tan(phi_r * m.pi/180)
    
    else:
        print('Gravel Size Sediment, Enter Sand or Smaller')
    
    return Tau_CR

#%%

#FERGUSON AND CHURCH (2006) PARTICLE SETTLING VELOCITY
#Input: diameter (m)
#Output: m/s | Input: micron
def SettlingVelocity(diameter):
    
    #Function Constants
    nu = 1*10**(-6)         #fluid kinematic vescosity (m2/s)
    g = 9.81                #gravitationa acceleration (m/s2)
    rho_f = 1.00            #fluid density (kg/m3)
    rho_s = 2.65            #solid density (kg/m3) 
    
    R = (rho_s - rho_f) / rho_f             #submerged grain specific gravity (unitless)
    diameter = diameter*(1*10**(-6))        #diameter converstion (m)

    #Stoke's Settling Velocity
    ws_stokes = (R*g*diameter**2)/(18*nu)
    Rep = (diameter*ws_stokes)/nu
    
    if Rep <= 1:
        ws = ws_stokes
  
    else:
        #Ferguson and Church (2006) Settling Velocity
        ws_fandc = (R*g*diameter**2)/(20*nu+(0.75*1.1*R*g*diameter**3)**0.5)
        Rep = (diameter*ws_fandc)/nu
        

        if 1 < Rep < 1000:
            ws = ws_fandc
            
        if Rep <= 1:
            ws = ws_stokes
            
        else:
            #Turbulent Flow Settling Velocity
            ws_turb = ((4*R*g*diameter)/(3*1.0))**(0.5)
            Rep = (diameter*ws_turb)/nu
            
            if 1000 < Rep:
                ws = ws_turb
            
            else:
                print ('Reynolds Particle Number (Rep) Error')
    
    return ws

#%%

#WIBERG & RUBIN (1989) BED LOAD CONCENTRATION
#Input: Bed Shear Stress (Pa); Critical Shear Stress (Pa)
#Output: m3/m3
def BedLoadConc(BedShearStress, CriticalShearStress):
    
    #Function Constants
    Con_MaxB = 0.65     #maximum bedload concentration (m3/m3) 
                   
    #Intermediate Calcualtions: T_Star (Sediment Mobility Fraction); E (Adjusted T_Star)
    T_star = BedShearStress / CriticalShearStress
    E = T_star - 1
    
    #Bed Load Concentration Calculation (m3/m3)
    BLConcentration = (0.045 * Con_MaxB * (E**0.93)) / (1 + 0.45*(E**0.93))
    
    return BLConcentration


#%%

#WIBERG & RUBIN (1989) SALTATION HEIGHT (BED LOAD HEIGHT)
#Input: diameter (microns); Bed Shear Stress (Pa); Critical Shear Stress (Pa)
#Output: m
def SaltationHeight(diameter, BedShearStress, CriticalShearStress):
    
    #Function Constants
    a1 = 0.68
    
    diameter = diameter*(1*10**(-6))    #diameter converstion (m)
    
    #Intermediate Calcualtions: T_Star (Sediment Mobility Fraction); a2 Fitted Equation
    T_star = BedShearStress / CriticalShearStress
    a2 = 0.0204 * (m.log(diameter))**2 + 0.0220*(m.log(diameter)) + 0.0709
    
    #Saltation Height Intermediate Calculation
    SaltH = diameter * ((a1 * T_star) / (1 + a2 * T_star))
    
    return SaltH

#%%

#VAN RIJN WAVE_RELATED BED SHEAR STRESS
#Input: depth (m)
#Output: Pa
def WaveShearStress(depth):   
  
    #USER DEFINED INPUTS
    L = 11                  #Wavelength (m)
    Hs = 0.2                #Significant Wave Heigh (m)
    Tp = 5                  #Peak Wave Period (sec)
    
    #Function Constants
    rho_f = 1000        #Density of water (kg/m3)
    nu = 10**(-6)       #Kintmatic viscosity of water (m2/s)
    
    #Intermediate Calculations: 
        #k (wave number); Uw (peak orbital velocity); Aw (Peak Orbital Excursion); fw (wave-related friction coefficient)
    
    k = (2 * m.pi) / L
    Uw = m.pi * Hs / (Tp * m.sinh(k*depth))
    Aw = (Tp / (2 * m.pi)) * Uw
    fw = 2 * ((Uw * Aw) / nu)**(-0.5)

    #Wave-Derived Bed Shear Stress Calculation
    tau_w = 0.25 * rho_f * fw * (Uw)**2

    return tau_w
