## The name of the program: FDEM_BfiPF
## Description

FDEM_BfiPF is a program that implements the numerical method for calculating the frequency-domain electromagnetic (EM) field of a grounded horizontal electric line (HEL) in an anisotropic horizontally layered medium.
In this method the desired field is sought by solving several axisymmetric problems using finite element approximations. 
Its advantages over known semi-analytical methods are manifested when 3D EM fields are calculated using the primary-secondary field approach and the simulated medium contains a sufficiently large number of 3D inhomogeneities. 
In the proposed method, the primary field of a HEL is represented as a sum of several axisymmetric fields whose origins are located at different points of the grounded electric line (field source). 
In this representation, the main components of the desired primary field are the fields of radial and point sources. 
An important characteristic of the proposed approach is that it completely eliminates the mutual subtraction of the components of the desired field and requires solving only two scalar axisymmetric problems.  

Program language: C++, Fortran  

---

##### The title of the manuscript: 
*Primary-secondary field approach to modeling frequency-domain 3D geoelectromagnetic field with FE approximation of primary field excited by a grounded horizontal electric line*  

##### The author details:

Yuri G. Soloveichik(a), Marina G. Persova(b*), Denis V. Vagin(c), and Yulia I. Koshkina(d)  

(a) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0003-2843-3214  
(b) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0003-1425-3538  
(c) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0002-2739-9563  
(d) Novosibirsk State Technical University, Novosibirsk, 630073, Russia, ORCID: 0000-0002-2684-4338  
* Corresponding author.  
E-mail address: mpersova@mail.ru (Marina G. Persova).  


## Requirements     
                    
**Features of the type of implementing computer or other computing device:**  a personal computer based on 64-bit Intel or AMD processors. 
 
**Operating system type and version:** Windows 10 and later. 
 
**Build (two ways):**
1. To compile using MinGW ("Minimalist GNU for Windows") for 64-bit versions of Windows (we use version 15.2.0), you must also install the Intel oneAPI Base Toolkit (we use version 2025.2), which includes the Intel oneAPI Math Kernel Library.

2. To build programs in Microsoft Visual Studio Community 2022 (we use version 17.14.18 (October 2025)) you must also install:
    - Intel Fortran Compiler (we use version 2025.3.0);
    - Intel oneAPI Base Toolkit (we use version 2025.2), which include the Intel oneAPI Math Kernel Library.

## How to build

Already built programs (modules) for Windows are in folder ['Folder_for_Calculations'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/Folder_for_Calculations) in https://github.com/Persova-LMDPHT/FDEM_BfiPF.   

If you want to build all programs (modules) yourself, download the contents of the ['CODE'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/CODE) folder from https://github.com/Persova-LMDPHT/FDEM_BfiPF. The ['CODE'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/CODE) folder contains the source code for all modules (each module is located in its own folder within the ['CODE'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/CODE) folder, with the same folder name as the module).   
To build all programs (modules) you can use one of the following methods:   

1. MinGW ("Minimalist GNU for Windows")  
 
   - To build modules using MinGW ("Minimalist GNU for Windows" for 64-bit versions of Windows (we use version 15.2.0)) you can use the following [makefile.bat](https://github.com/Persova-LMDPHT/FDEM_BfiPF/blob/731ce1cca237ad1e70d51f6a0f91fb799bf2d7e6/CODE/makefile.bat) as an instruction (you must first install Intel oneAPI (we use version 2025.2), which should include the Intel oneAPI Math Kernel Library).

   - In [makefile.bat](https://github.com/Persova-LMDPHT/FDEM_BfiPF/blob/731ce1cca237ad1e70d51f6a0f91fb799bf2d7e6/CODE/makefile.bat), first change the path to the folder with *.h MKL files (variable %PATH_TO_MKL_INCLUDES%) and the path to the folder with MKL libraries (variable %PATH_TO_MKL_LIBS%) to the paths corresponding to your computer.

   - Go to the folder containing the source code for all modules (each module is located in its own folder, the folder name matches the module name) and run [makefile.bat](https://github.com/Persova-LMDPHT/FDEM_BfiPF/blob/731ce1cca237ad1e70d51f6a0f91fb799bf2d7e6/CODE/makefile.bat) there. The compiler messages are displayed in !logg++.txt or !loggfortran.txt files in the folder containing the corresponding module's source code.

   - The [makefile.bat](https://github.com/Persova-LMDPHT/FDEM_BfiPF/blob/731ce1cca237ad1e70d51f6a0f91fb799bf2d7e6/CODE/makefile.bat) commands will create the folder 'Folder_for_Calculations' containing the CalcStarter.exe executable file and the subfolder 'Folder_for_Calculations\Modules' containing all 18 necessary executable files. When running programs (modules), DLL files may be required, which must be placed to the 'Modules' folder. For user convenience, we have placed the main DLL files in the ['Folder_for_Calculations/Modules/'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/Folder_for_Calculations/Modules) folder in https://github.com/Persova-LMDPHT/FDEM_BfiPF. 

If everything is done correctly, your 'Folder_for_Calculations' is ready to use!   

To perform the calculation, add input files (they are described below) from ['CalculationExample'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/89293a537e8dbebc5c7704a6e67aa143482e42f4/CalculationExample) to the 'Folder_for_Calculations' and run CalcStarter.exe.

2. Microsoft Visual Studio   

To build programs (modules) in Microsoft Visual Studio, please see section ‘Requirements’ and read the file ['Instruction_for_Modules_Build_in_Visual_Studio.pdf'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/blob/aecd255a0692bd37c58a5bb66d2b673e3c49d496/Instruction_for_Modules_Build_in_Visual%20Studio.pdf).


## How to run calculation

To run calculation follow the steps **(do not use special characters and too long paths to the folder for calculations)**:
1. Build all programs from folder ['CODE'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/CODE) and go to step 2 
   or use a ready-made folder for calculation  with already built programs (folder ['Folder for Calculations'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/Folder_for_Calculations)) and go to step 7.
2. Create empty folder for calculation running (let's denote it as 'Test' folder). *Do not use special characters and too long paths to the folder for calculations!*
3. Create 'Modules' folder inside of 'Test' folder.
4. Put all calculation programs (av.exe, bound.exe, CalcFreq.exe, CalcHarm2D_Ax.exe, CalcHarm2DHEL.exe, CalcHarm2DHEL_AV.exe, CalcHarm2DHEL_U.exe, CalcHarm3D.exe, Mesh2D_FD.exe, OutputSmooth2DHarm_Ax.exe, OutputSmooth2DHarm_U.exe, OutputSmoothAV2DHarm_Ax.exe, OutputSmoothAV2DHarm_Er.exe, OutputSmoothAV2DHarm_Ez.exe, RegularMeshBuilder.exe, SumHarm2D3D.exe, u.exe, UnloadAnomalHarm.exe) to 'Modules'
5. (If needed) Put required dll files to 'Modules' folder.
6. Put CalcStarter.exe to 'Test' folder.
*The example of the ready-made folder for calculation with already built programs is in folder ['Folder for Calculations'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/aecd255a0692bd37c58a5bb66d2b673e3c49d496/Folder_for_Calculations) in https://github.com/Persova-LMDPHT/FDEM_BfiPF*
7. Put all prepared input text files (gen, layers, objects, recE, recB, settings.cfg) to 'Test' folder.   You can take already prepared input files from ['CalculationExample'](https://github.com/Persova-LMDPHT/FDEM_BfiPF/tree/89293a537e8dbebc5c7704a6e67aa143482e42f4/CalculationExample).
8. Run CalcStarter.exe and wait for it to complete.
9. Look for result file (e2d.1, e3d.1) in 'Results' folder inside of 'Test' folder.

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