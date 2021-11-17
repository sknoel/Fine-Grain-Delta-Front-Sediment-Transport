
"""
@author: sknoel
"""

#PYTHON IMPORTS
from pathlib import Path
from matplotlib import pyplot as plt

import math as m
import pandas as pd
import os

#Local Function Directory
Directory = os.chdir("INSERT FILE PATH HERE \Functions")
from Functions import Tau_DepthSlope, Tau_LotW, Tau_CR, SettlingVelocity, BedLoadConc, SaltationHeight, WaveShearStress

#%%
#MODIFICATIONS
#Tidal Increase (m)
Tides = True
TideDepth = 0.3

#Wave Influence
Waves = False

print ('Model Parameters:', '(1) Tides =', Tides, '(1.1) Tide Depth (m) =', TideDepth, '(2) Waves =', Waves)

#%%
#READ IN: User Input from CSV
model_folder =  Path("INSERT FILE PATH HERE \Model Input Files")
BinnedUserInput = model_folder/"Rouse Profile Bins.csv"
SingularUserInput = model_folder/"Rouse Profile Inputs.csv"

#Binned User Input to lists
BUI = pd.read_csv(BinnedUserInput)
BUI = BUI.dropna(axis=1, how='all')
BUI = BUI.dropna(axis=0, how='all')

Bin_Length = BUI.iloc[0].values.tolist()[1:]
Depth = BUI.iloc[1].values.tolist()[1:]
Velocity = BUI.iloc[2].values.tolist()[1:]

if Waves:
    Wave_Attenuation = BUI.iloc[3].values.tolist()[1:]
else:
    Wave_Attenuation = [90]*len(Bin_Length)

if Tides:
    i = 0
    print('Tides Used in Calculation. Tide value input: ', TideDepth)
    
    while i < (len(Bin_Length)):
        Depth[i] = Depth[i] + TideDepth
        i = i + 1
        continue

#Singular User Inputs to lists
SUI = pd.read_csv(SingularUserInput, header=None)
SUI = SUI.dropna(axis=1, how='all')
SUI = SUI.dropna(axis=0, how='all')

Grain_Sizes = SUI.iloc[0].values.tolist()[1:]
N = len(Grain_Sizes)

#%%
#Reach Slope Calculation for Bed Shear Stress for first itteration
ReachSlope = []
i = 0
while i < (len(Depth) - 1):
    ReachSlope.append((Depth[i+1] - Depth[i])/Bin_Length[i])
    i = i + 1
    continue

ReachSlope.append(ReachSlope[i-1])

#%%
#BEGIN MODEL LOOP CALCULATIONS
    #for loop = grain size inputs
    #first while loop = bins
    #second while loop = itterative bed shear stress
    #Thirs while loop = numerical method for partical differenctional equation

#CONSTANTS
k = 0.4             #von Karman constant (1930), unitless
rho_f = 1000        #Fluid density (kg/m3)
    
x = N
y = len(Bin_Length)
Rouse_Profile = []
UCH = []

j = 0
for GS in Grain_Sizes:
    Rouse_Profile.append([])
    UCH.append([])
    
    ws = SettlingVelocity(Grain_Sizes[j])
    TauCR = Tau_CR(Grain_Sizes[j])
    u_star = ( TauCR / rho_f )**0.5
    
    RouseNumber = ws / (k * u_star)
    #print('RouseNumber', RouseNumber)

    i = 0
    while i < len(Bin_Length):
        
        #Itterative Bed Shear Stress  
        Tau_Bed = Tau_DepthSlope(Depth[i], ReachSlope[i])

        while True:
            
            C_Ref = BedLoadConc(Tau_Bed, TauCR)
            SaltH = SaltationHeight(Grain_Sizes[j], Tau_Bed, TauCR)
            TB = Tau_Bed
            Tau_Bed = Tau_LotW(Velocity[i], Depth[i], SaltH)
            
            if abs( Tau_Bed - TB ) / Tau_Bed < 0.001 or m.isnan(Tau_Bed):   #Tolerance Value
                break        
            
            #Final Bed Shear Stress Value
            Tau_Bed = Tau_Bed + WaveShearStress(Depth[i]) * (m.cos( Wave_Attenuation[i] * m.pi/180)) 
            continue
        
        C_Grain = (1/N) * C_Ref * (Depth[i]**(-1*RouseNumber) * (SaltH * (Depth[i] ** RouseNumber) - Depth[i] * (SaltH ** RouseNumber))) / ((RouseNumber - 1) * (Depth[i] - SaltH))
        
        top = Velocity[i] * C_Grain * Depth[i]              
        UCH[j].append(top)
            
        i = i + 1
        continue
   
    i = 0
    while i < len(Bin_Length):
        
        if i == 0:
            Rouse_Profile[j].append(( UCH[j][i + 1] - UCH[j][i] ) / ( Bin_Length[i + 1] ))
    
        elif i == ( len(Bin_Length) - 1 ):
            Rouse_Profile[j].append(( UCH[j][i] - UCH[j][i - 1] ) / ( Bin_Length[i] ))
        
        else:
            Rouse_Profile[j].append(( UCH[j][i + 1] - UCH[j][i - 1] ) / ( Bin_Length[i] + Bin_Length[i + 1] ))
        
        i = i + 1
        continue

    j = j + 1

#%%
#GENERATED GRAPHICS
#Produces single figure with all (maximum 2) mass conservation profiles.

#Cummulative distance for graphcs generation
Cumm_Bin=[]
Cumm_Bin.append(Bin_Length[0])
i = 1
while i < (len(Bin_Length)):
    Cumm_Bin.append(Cumm_Bin[i-1]+ Bin_Length[i])
    i = i + 1

SILT = Rouse_Profile[0]
SAND = Rouse_Profile[1]

ZERO = []
i = 0
while i < (len(Bin_Length)):
    ZERO.append(0)
    i = i + 1

figure, axis = plt.subplots(2,1)
plt.subplots_adjust(hspace = 0.5)

figure.suptitle('Rouse Concentration Profile Based Model', fontsize=14)
figure.text(0.5, 0.04, 'Distance from Subaqueous Channel, x (m)', ha='center')
figure.text(0.04, 0.5, 'Mass Conservation Partial Derivative, Ri(x)', va='center', rotation='vertical')

axis[0].plot(Cumm_Bin, SILT, color='chocolate', linewidth=4)
axis[0].plot(Cumm_Bin, ZERO, color='black', linewidth=0.5)
axis[0].set_title('Silt')
axis[0].grid(which='major', alpha=0.6)
axis[0].grid(which='minor', alpha=0.6)
axis[0].set_xlim([203, max(Cumm_Bin)])

axis[1].plot(Cumm_Bin, SAND, color='goldenrod', linewidth=4)
axis[1].plot(Cumm_Bin, ZERO, color='black', linewidth=0.5)
axis[1].set_title('Sand')
axis[1].grid(which='major', alpha=0.6)
axis[1].set_xlim([203, max(Cumm_Bin)])
