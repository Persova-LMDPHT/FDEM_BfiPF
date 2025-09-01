The name of the program: FDEM_BfiPF
---------------------------------
FDEM_BfiPF is a program that implements the numerical method for calculating the frequency-domain electromagnetic (EM) field of a grounded horizontal electric line (HEL) in an anisotropic horizontally layered medium.
In this method the desired field is sought by solving several axisymmetric problems using finite element approximations. 
Its advantages over known semi-analytical methods are manifested when 3D EM fields are calculated using the primary-secondary field approach and the simulated medium contains a sufficiently large number of 3D inhomogeneities. 
In the proposed method, the primary field of a HEL is represented as a sum of several axisymmetric fields whose origins are located at different points of the grounded electric line (field source). 
In this representation, the main components of the desired primary field are the fields of radial and point sources. 
An important characteristic of the proposed approach is that it completely eliminates the mutual subtraction of the components of the desired field and requires solving only two scalar axisymmetric problems. 
---------------------------------
The title of the manuscript:
Primary-secondary field approach to modeling frequency-domain 3D geoelectromagnetic field with FE approximation of primary field excited by a grounded horizontal electric line
--------------------------------- 
The author details:
Yuri G. Soloveichik(a), Marina G. Persova(b*), Denis V. Vagin(c), and Yulia I. Koshkina(d)

(a) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0003-2843-3214
(b) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0003-1425-3538
(c) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0002-2739-9563
(d) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0002-2684-4338
* Corresponding author.
E-mail address: mpersova@mail.ru (Marina G. Persova).


//======== How to run calculation ===================================											  
To run calculation follow the steps (do not use long paths to the directory for calculation!!!):
1. Build all programs from folder 'CODE' and go to item 2 
   or use a ready-made folder for calculation  with already assembled programs (folder 'Folder for Calculations') and go to item 7.
2. Create empty directory for calculation running (let's denote it as Home directory).
3. Create 'Modules' directory inside of Home directory.
4. Put all calculation programs (av.exe, bound.exe, CalcFreq.exe, CalcHarm2D_Ax.exe, CalcHarm2DHEL.exe, CalcHarm2DHEL_AV.exe, CalcHarm2DHEL_U.exe, CalcHarm3D.exe, Mesh2D_FD.exe, OutputSmooth2DHarm_Ax.exe, OutputSmooth2DHarm_U.exe, OutputSmoothAV2DHarm_Ax.exe, OutputSmoothAV2DHarm_Er.exe, OutputSmoothAV2DHarm_Ez.exe, RegularMeshBuilder.exe, SumHarm2D3D.exe, u.exe, UnloadAnomalHarm.exe) to 'Modules'
5. (If needed) Put required dll files to ' Modules' directory.
6. Put CalcStarter.exe to Home directory.
7. Put all prepared text files (gen, layers, objects, recE, recB, settings.cfg) to Home directory.
8. Run CalcStarter.exe and wait for it to complete.
9. Look for result file (e2d.1, e3d.1) in 'Results' directory inside of Home directory.

All dimensions are in meters. In the calculations, the transmitter current is taken equal to 1 A.



//=========Input files=======================

1. layers
This file contains the resistivity of layers defined in the model. Do not include air in this file. The air layer will be added automatically

Comment.

Note that when adding an air layer to the computational domain, we use sigma=1e-8 S/m. This is due to the fact that in surface integrals, 
sigma is in the denominator (therefore, sigma=0 is unacceptable). At the same time, studies show (see Fig. 5 and the description to it) 
that it is better to take the contribution to the matrix from the surface integral on the air side (which is ensured by the weights calculated using formula (111)). 
In this case, increasing the sigma value, for example, to 1e-6 S/m, can lead to an increase in error due to the fact that it does not correspond to air, 
and decreasing the sigma value, for example, to 1e-12 S/m, also leads to an increase in error and requires decreasing the mesh step hz to 0.001 m. 
Therefore, it is recommended to use sigma=1e-8 S/m and hz=0.1 m. The correctness of this choice is confirmed by comparison with the Dipole1D code.
              

File format:
(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
<Layers count>
For each layer:
<Upper bound Z,m> <Conductivity_l, S/m> <Conductivity_z, S/m>

=============================================================================================
2. objects
This file contains 3D objects defined as parallelepipeds.
(!) Note that the Z axis is directed upwards (see Fig. 6a of the article). 

File format:
<Objects count>
For each polygon:
<X0, m> <X1, m> <Y0, m> <Y1, m> <Z0, m> <Z1, m> <Conductivity_l, S/m> <Conductivity_z, S/m>

=============================================================================================
3. gen
This file contains transmitter (AB) line placement.
(!) Note that the Z axis is directed upwards (see Fig. 6a of the article). 

File format:
<AB_Z, m>
<A_X, m> <A_Y, m> <B_X, m> <B_Y, m>

=============================================================================================
4. recE
This file contains common list of receiver placements.
(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).

File format:
<Receivers count>
For each receiver:
<X, m> <Y, m> <Z, m>

=============================================================================================
5. recB                                                                                               
This file contains common list of receiver placements.
(!) Note that the Z axis is directed upwards (see Fig. 6a of the article).
                                                 
                                                                                                      
File format:                                                                                          
<Receivers count>                                                                                     
For each receiver:                                                                                    
<X, m> <Y, m> <Z, m>                                                                                           
                                                                                                      
=============================================================================================         
6. settings.cfg
This file contains mesh settings and some calculation settings

File format:

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
                                                       
Method for 2D (0 - Bfi-formulation, 1 - AV-formulation)             
                                                       
SLAE solver for 3D problem (0 - direct solver, 1 - iterative solver)     
<Minimum residual value for exiting the iterative process (for 3D problem if the "iterative solver" is used)>
<Maximum allowed number of iterations for exiting the iterative process (for 3D problem if the "iterative solver" is used)>
                                                       


//============ Results files ===================================
'e2d.1', 'e3d.1' files contain the signals in all receivers in file 'recE'.
<Re_Ex, B/m> <Re_Ey, B/m>  <Re_Ez, B/m> <Im_Ex, B/m> <Im_Ey, B/m> <Im_Ez, B/m>

'b2d.1', 'b3d.1' files contain the signals in all receivers in file 'recB'.  
<Re_Bx, T> <Re_By, T>  <Re_Bz, T> <Im_Bx, T> <Im_By, T> <Im_Bz, T>                             

