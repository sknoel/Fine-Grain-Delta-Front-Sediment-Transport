
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
from Functions import Tau_DepthSlope, Tau_LotW, Tau_CR, SettlingVelocity, SaltationHeight, WaveShearStress

#%%
#USER DEFINED INPUTS
r0 = 0.25          #bed roughness height (unitless)

#MODIFICATIONS
#Tidal Increase (m)
Tides = False
TideDepth = 0.3

#Wave Influence
Waves = False

print ('Model Parameters:', '(1) Tides =', Tides, '(1.1) Tide Depth (m) =', TideDepth, '(2) Waves =', Waves)

#%%
#READ IN: User Input from CSV
model_folder =  Path("INSERT FILE PATH HERE \Model Input Files")
BinnedUserInput = model_folder/"Advection Settling Bins.csv"
SingularUserInput = model_folder/"Advection Settling Inputs.csv"

#Binned User Input to lists
BUI = pd.read_csv(BinnedUserInput, header=None)
BUI = BUI.dropna(axis=1, how='all')
BUI = BUI.dropna(axis=0, how='all')

Bin_No = BUI.iloc[0].values.tolist()[1:]
Bin_Length = BUI.iloc[1].values.tolist()[1:]
Depth = BUI.iloc[2].values.tolist()[1:]
Velocity = BUI.iloc[3].values.tolist()[1:]

if Waves:
    Wave_Attenuation = BUI.iloc[4].values.tolist()[1:]
else:
    Wave_Attenuation = [90]*len(Bin_No)

#Singular User Inputs to lists
SUI = pd.read_csv(SingularUserInput, header=None)
SUI = SUI.dropna(axis=1, how='all')
SUI = SUI.dropna(axis=0, how='all')

Grain_Sizes = SUI.iloc[0].values.tolist()[1:]
Starting_Depth = SUI.iloc[1].values.tolist()[1:]

#%%
#Reach Slope Calculation for Bed Shear Stress for first itteration
ReachSlope = []
i = 0
while i < (len(Depth) - 1):
    ReachSlope.append((Depth[i+1] - Depth[i])/Bin_Length[i])
    i = i + 1
    continue
ReachSlope.append(ReachSlope[i-1])

if Tides:
    i = 0
    print('Tides Used in Calculation. Tide value input: ', TideDepth)
    
    while i < (len(Bin_Length)):
        Depth[i] = Depth[i] + TideDepth
        i = i + 1
        continue

#%%
#BEGIN MODEL LOOP CALCULATIONS
    #for loop = grain size inputs
    #first while loop = bins
    #second while loop = itterative bed shear stress
    
Particle_Depth = []
j = 0
for GS in Grain_Sizes:
    Particle_Depth.append([])    
    Particle_Depth[j].append(Starting_Depth[j])
    
    ws = SettlingVelocity(GS)
    Tau_CR = Tau_CR(Grain_Sizes[j])
    
    i = 0
    while i < len(Bin_Length):
 
        if (not Particle_Depth[j][i] == Depth[i-1]) :
            z = ( (Bin_Length[i] * r0 * ws) / Velocity[i] ) + Particle_Depth[j][i]          
        
        else: 
            z = Depth[i]
    
        if z > Depth[i]:
            z = Depth[i]
            
            #Itterative Bed Shear Stress           
            Tau_Bed = Tau_DepthSlope(Depth[i], ReachSlope[i])

            while True:

                SaltH = SaltationHeight(Grain_Sizes[j], Tau_Bed, Tau_CR)
                TB = Tau_Bed
                Tau_Bed = Tau_LotW(Velocity[i], Depth[i], SaltH)
                
                if abs( Tau_Bed - TB ) / Tau_Bed < 0.001 :   #Tolerance Value 
                    break        

            #Final Bed Shear Stress Value
            Tau_Bed = Tau_Bed + WaveShearStress(Depth[i]) * (m.cos( Wave_Attenuation[i] * m.pi/180))

            if Tau_Bed < Tau_CR:
                break
    
        Particle_Depth[j].append(z)
        i = i + 1
        continue

    j = j + 1

#%%
#GENERATED GRAPHICS
#Produces single figure with all (maximum 7 distinct) grain sizes.
#Maximum 7 distinct grain sizes in excisting graphics format (additional can be added as "elif" statements)

#Cummulative distance for graphic generation
Cumm_Bin = []
Cumm_Bin.append(Bin_Length[0])
i = 1
while i<( len(Bin_Length) ):
    Cumm_Bin.append(Cumm_Bin[i-1] + Bin_Length[i])
    i = i + 1

#Grain size depth plot. Grain size dependent distance traveled.
plt.figure(1)
plt.title('Advection Settling Model')
plt.ylabel('Depth, h (m)')
plt.xlabel('Distance from Subaqueous Channel, x (m)')
plt.ylim(max(Depth), 0)
plt.xlim(0, max(Cumm_Bin))
plt.plot(Cumm_Bin, Depth, color='wheat', linestyle='-', linewidth=8,  label='Bed Level')

j = 0
for GS in Grain_Sizes:
    Plot_Bins = [0]
    
    i = 0
    while i < len(Particle_Depth[j]) - 1 :
        Plot_Bins.append(Bin_Length[i]+Plot_Bins[i])
        i = i + 1
        continue

    if j == 0:
        plt.plot(Plot_Bins, Particle_Depth[j], color='red', linestyle='dashdot', linewidth=2, label=Grain_Sizes[j])

    elif j == 1:
        plt.plot(Plot_Bins, Particle_Depth[j], color='blue', linestyle='dashed', linewidth=2, label=Grain_Sizes[j])
    
    elif j == 2:
        plt.plot(Plot_Bins, Particle_Depth[j], color='indigo', linestyle='dotted', linewidth=2, label=Grain_Sizes[j])
        
    elif j == 3:
        plt.plot(Plot_Bins, Particle_Depth[j], color='green', linestyle=(0, (5, 5)), linewidth=2, label=Grain_Sizes[j])
    
    elif j == 4:
        plt.plot(Plot_Bins, Particle_Depth[j], color='grey', linestyle=(0, (1,1)), linewidth=2, label=Grain_Sizes[j])
        
    elif j ==5:
        plt.plot(Plot_Bins, Particle_Depth[j], color='orangered', linestyle=(0, (3, 1, 1, 1)), linewidth=2, label=Grain_Sizes[j])
     
    elif j > 5:
        plt.plot(Plot_Bins, Particle_Depth[j], color='yellow', linestyle=(0, (5, 1)), linewidth=2, label=Grain_Sizes[j])

    plt.legend(loc="best", prop={'size':8})
    j = j + 1

