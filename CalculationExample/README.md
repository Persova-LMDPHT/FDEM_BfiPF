## This folder contains the folders with examples of input files and calculation results for various geoelectric models.

Descriptions of the input and result files are given below.

## Input files (layers, objects, gen, recE, recB, settings.cfg)
##### ♦ All lengths are in meters. In the calculations, the transmitter current is taken equal to 1 A.
#### 1. layers
This file contains the resistivity of layers defined in the model. **Do not include air in this file. The air layer will be added automatically**

##### ♦ Comment:

*Note that when adding an air layer to the computational domain, we use sigma=1e-8 S/m. This is due to the fact that in surface integrals, 
sigma is in the denominator (therefore, sigma=0 is unacceptable). At the same time, studies show (see Fig. 5 and the description to it) 
that it is better to take the contribution to the matrix from the surface integral on the air side (which is ensured by the weights calculated using formula (111)). 
In this case, increasing the sigma value, for example, to 1e-6 S/m, can lead to an increase in error due to the fact that it does not correspond to air, 
and decreasing the sigma value, for example, to 1e-12 S/m, also leads to an increase in error and requires decreasing the mesh step hz to 0.001 m. 
Therefore, it is recommended to use sigma=1e-8 S/m and hz=0.1 m. The correctness of this choice is confirmed by comparison with the Dipole1D code.*
              

File format:
```(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<Layers count>
For each layer:
<Upper bound Z,m> <Conductivity_l, S/m> <Conductivity_z, S/m>
```

#### 2. objects
This file contains 3D objects defined as parallelepipeds.
 
File format:
```(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<Objects count>
For each polygon:
<X0, m> <X1, m> <Y0, m> <Y1, m> <Z0, m> <Z1, m> <Conductivity_l, S/m> <Conductivity_z, S/m>
```
#### 3. gen
This file contains transmitter (AB) line placement.

File format:
```(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<AB_Z, m>
<A_X, m> <A_Y, m> <B_X, m> <B_Y, m>
```

#### 4. recE
This file contains common list of receiver placements.

File format:
```(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<Receivers count>
For each receiver:
<X, m> <Y, m> <Z, m>
```
#### 5. recB 
This file contains common list of receiver placements.

File format:
```(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<Receivers count> 
For each receiver:
<X, m> <Y, m> <Z, m> 
```

#### 6. settings.cfg
This file contains mesh settings and some calculation settings

File format:
```
Mesh 3D section:
<Initial step X, Y> <Initial step Z>
<Sparce coefficient X, Y> <Sparce coefficient Z>
<Gap between receivers bounding box and sparced mesh area X, Y> <Gap between receivers bounding box and sparced mesh area Z>
<Distance from receivers bounding box+gap to the boundary of the calculation domain X, Y> <Distance from receivers bounding box+gap to the boundary of the calculation domain Z>
<Number of refinements around recevers>

<Approximation dipoles count for transmitter line>  
Mesh 2D section:
<The distance to remote boundary>
For Bfi-formulation: <minimum step along the r-axis> <minimum step along the z-axis> <step increase coefficient along the r-axis> <step increase coefficient along the z-axis inside layers> <step increase coefficient along the z-axis when moving to the boundaries of the domain>
For  AV-formulation: <minimum step along the r-axis> <minimum step along the z-axis> <step increase coefficient along the r-axis> <step increase coefficient along the z-axis inside layers> <step increase coefficient along the z-axis when moving to the boundaries of the domain>

<Frequency>

Calculations section:
                                                       
<Threads count for parallel calculation>
                                                       
Method for 2D:<0> - Bfi-formulation, <1> - AV-formulation
                                                       
SLAE solver for 3D problem: <0> - direct solver, <1> - iterative solver     
<Minimum residual value for exiting the iterative process (for 3D problem if the "iterative solver" is used)>
<Maximum allowed number of iterations for exiting the iterative process (for 3D problem if the "iterative solver" is used)>
```

## Results files 
**'e2d.1'**, **'e3d.1'** files contain the signals in all receivers in file 'recE'.
File format:
```
<Re_Ex, B/m> <Re_Ey, B/m>  <Re_Ez, B/m> <Im_Ex, B/m> <Im_Ey, B/m> <Im_Ez, B/m>
```

**'b2d.1'**, **'b3d.1'** files contain the signals in all receivers in file 'recB'. 
File format: 
```
<Re_Bx, T> <Re_By, T>  <Re_Bz, T> <Im_Bx, T> <Im_By, T> <Im_Bz, T>
```
The result files are in the 'Results' folder inside the folder for calculation.