# University of New Orleans
# Department of Earth and Environmental Sciences

FINE GRAINED DELTA FRONT SEDIMENT TRANSPORT
M.S. Earth and Environmental Science
Sarah Noel, December 2021

Description of Models:
The enclosed code contains two distinct models (1) an augmented Advection Settling and (2) Mass Balance Rouse Profile. 
Referenced functions are catalogued in the funnction.py file with commented method referenced, inputs, and outputs. 

(1) Advection Settling Model
The augmented advection settling model follows the path of specific grain size as it falls though the water column.
Once the grain reaches the bed, the bed shear stress conditions are evaluated to determine if the grain can be carried
as bed load or will be deposited. This depositional check augments the interpretation of advected sediment transport.

(2) Mass Balance Rouse Profile Model
The mass balance Rouse profile model considers the fluid stress conditions, grain size, and water depth to 
determine the erosional or depositional state for a specific grain size fraction at a given location. A negative value
indicates a depositional state, while a positive value indicated an erosional state. 


Installation and Use Guide:
1) Unzip and place nested folders into a stable location.
2) Update the following file paths:
        (1) Advection Settling lines #15, 34
        (2) Mass Balance Rouse Profile line #15, 31