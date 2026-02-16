md Folder_for_Calculations
md Folder_for_Calculations\Modules 
set "PATHTOCALCFOLDER=..\Folder_for_Calculations\"

set "PATH_TO_MKL_INCLUDES=C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2/include"
set "PATH_TO_MKL_LIBS=C:/Program Files (x86)/Intel/oneAPI/mkl/2025.2/lib"
:: build Module av.exe 
cd ./av
"g++.exe" -static -O2 -Wall -std=c++14 -I "%PATH_TO_MKL_INCLUDES%" -L "%PATH_TO_MKL_LIBS%" bound_cond_2d.cpp FormatConverter.cpp For_Solvers.cpp global_slae_2d.cpp in_out.cpp Local_matrix_rz.cpp MeshRZ.cpp pardiso.cpp Portret.cpp stdafx.cpp VEL2D.cpp VelHarm2d.cpp -l mkl_intel_lp64 -l mkl_intel_lp64_dll -l mkl_core_dll -l mkl_intel_thread_dll -l mkl_rt -l mkl_sequential_dll -o av 2>!logg++.txt 
copy av.exe  %PATHTOCALCFOLDER%Modules 
cd ../
:: build Module bound.exe 
cd ./bound
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp stdafx.cpp task2d.cpp -o bound 2>!logg++.txt      
copy bound.exe  %PATHTOCALCFOLDER%Modules  
cd ../ 
:: build Module CalcFreq.exe
cd ./CalcFreq
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp -o CalcFreq 2>!logg++.txt 
copy CalcFreq.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module CalcHarm2DHEL.exe
cd ./CalcHarm2DHEL 
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp -o CalcHarm2DHEL 2>!logg++.txt
copy CalcHarm2DHEL.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module CalcHarm2DHEL_AV.exe
cd ./CalcHarm2DHEL_AV       
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp -o CalcHarm2DHEL_AV 2>!logg++.txt
copy CalcHarm2DHEL_AV.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module CalcHarm2DHEL_U.exe
cd ./CalcHarm2DHEL_U        
"g++.exe" -static -O1 -Wall -std=c++14  main.cpp -o CalcHarm2DHEL_U 2>!logg++.txt
copy CalcHarm2DHEL_U.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module CalcHarm2D_Ax.exe 
cd ./CalcHarm2D_Ax 
"g++.exe" -static -O2 -Wall -std=c++14 -I "%PATH_TO_MKL_INCLUDES%" -L "%PATH_TO_MKL_LIBS%" bound_cond_2d.cpp FormatConverter.cpp For_Solvers.cpp global_slae_2d.cpp Harm2D.cpp Harm2dLine.cpp in_out.cpp iobinary.cpp localcontributions.cpp pardiso.cpp Portret.cpp rect_local_matrix_2d.cpp rect_postproc_2d.cpp stdafx.cpp -l mkl_intel_lp64 -l mkl_intel_lp64_dll -l mkl_core_dll -l mkl_intel_thread_dll -l mkl_rt -l mkl_sequential_dll -o CalcHarm2D_Ax 2> !logg++.txt 
copy CalcHarm2D_Ax.exe  %PATHTOCALCFOLDER%Modules 
cd ../
:: build Module CalcHarm3D.exe  
cd ./CalcHarm3D
"g++.exe" -static -O2 -Wall -std=c++14 -fopenmp -I "%PATH_TO_MKL_INCLUDES%" -L "%PATH_TO_MKL_LIBS%" base_solver.cpp block_2x2_solver.cpp bound_cond_vec_harm.cpp CheckInHex.cpp ControlOMP.cpp ElemNeib.cpp FoldedPreconditioner.cpp FormatConverter.cpp For_Solvers.cpp give_out_vec_mt.cpp HarmVect3d.cpp Hex_Local_Matrix.cpp in_out.cpp iobinary.cpp MRS.cpp OutputArbitrary.cpp OutputResultant3d.cpp pardiso.cpp pcocr.cpp pcocr_rci.cpp Portret.cpp rci.cpp stdafx.cpp Subdomain.cpp T_Brick.cpp t_global_slae.cpp T_Mapping.cpp T_Portrait.cpp vec_prep_data.cpp -lpthread -l mkl_intel_lp64 -l mkl_intel_lp64_dll -l mkl_core_dll -l mkl_intel_thread_dll -l mkl_rt -l mkl_sequential_dll  -o CalcHarm3D  2> !logg++.txt
copy CalcHarm3D.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module Mesh2D_FD.exe  
cd ./Mesh2D_FD
"gfortran.exe" -static -O2 -Wall inputF.for Convert_data.for Meah1DZ.for Model_Correct.for Make_NormArea.for -o Mesh2D_FD 2> !loggfortran.txt 
copy Mesh2D_FD.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module OutputSmooth2DHarm_Ax.exe  
cd ./OutputSmooth2DHarm_Ax 
"g++.exe" -static -O2 -Wall -std=c++14 -fopenmp SimpleOutput2D.cpp SmoothOutput2D.cpp stdafx.cpp task2d.cpp main.cpp -lpthread  -o OutputSmooth2DHarm_Ax  2> !logg++.txt 
copy OutputSmooth2DHarm_Ax.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module OutputSmooth2DHarm_U.exe  
cd ./OutputSmooth2DHarm_U 
"g++.exe" -static -O2 -Wall -std=c++14 -fopenmp main.cpp SmoothOutput2D.cpp stdafx.cpp task2d.cpp -lpthread -o OutputSmooth2DHarm_U  2> !logg++.txt  
copy OutputSmooth2DHarm_U.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module OutputSmoothAV2DHarm_Ax.exe  
cd ./OutputSmoothAV2DHarm_Ax
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp SimpleOutput2D.cpp SmoothOutput2D.cpp stdafx.cpp task2d.cpp -o OutputSmoothAV2DHarm_Ax  2> !logg++.txt
copy OutputSmoothAV2DHarm_Ax.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module OutputSmoothAV2DHarm_Er.exe  
cd ./OutputSmoothAV2DHarm_Er
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp SimpleOutput2D.cpp SmoothOutput2D.cpp stdafx.cpp task2d.cpp -o OutputSmoothAV2DHarm_Er  2> !logg++.txt
copy OutputSmoothAV2DHarm_Er.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module OutputSmoothAV2DHarm_Ez.exe  
cd ./OutputSmoothAV2DHarm_Ez
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp SimpleOutput2D.cpp SmoothOutput2D.cpp stdafx.cpp task2d.cpp -o OutputSmoothAV2DHarm_Ez  2> !logg++.txt
copy OutputSmoothAV2DHarm_Ez.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module RegularMeshBuilder.exe  
cd ./RegularMeshBuilder
"g++.exe" -static -O2 -Wall -std=c++14 BasicOperations.cpp Logging.cpp main.cpp Mesh1D.cpp Mesh2D.cpp Mesh3D.cpp Model.cpp Processing.cpp ReadWrite.cpp -o RegularMeshBuilder  2> !logg++.txt
copy RegularMeshBuilder.exe  %PATHTOCALCFOLDER%Modules  
cd ../
:: build Module SumHarm2D3D.exe  
cd ./SumHarm2D3D
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp -o SumHarm2D3D 2>!logg++.txt
copy SumHarm2D3D.exe  %PATHTOCALCFOLDER%Modules
cd ../ 
:: build Module u.exe  
cd ./u
"g++.exe" -static -O2 -Wall -std=c++14 -I "%PATH_TO_MKL_INCLUDES%" -L "%PATH_TO_MKL_LIBS%" bound_cond_2d.cpp FormatConverter.cpp For_Solvers.cpp global_slae_2d.cpp in_out.cpp Local_matrix_rz.cpp MeshRZ.cpp pardiso.cpp Portret.cpp stdafx.cpp VEL2D.cpp VelHarm2d.cpp -l mkl_intel_lp64 -l mkl_intel_lp64_dll -l mkl_core_dll -l mkl_intel_thread_dll -l mkl_rt -l mkl_sequential_dll -o u  2> !logg++.txt
copy u.exe  %PATHTOCALCFOLDER%Modules 
cd ../
:: build Module UnloadAnomalHarm.exe  
cd ./UnloadAnomalHarm 
"g++.exe" -static -O2 -Wall -std=c++14 GeoCalculatorSP.cpp in_out.cpp stdafx.cpp Tensor.cpp time_approx.cpp T_Mapping.cpp vec_prep_data.cpp -o UnloadAnomalHarm  2> !logg++.txt
copy UnloadAnomalHarm.exe  %PATHTOCALCFOLDER%Modules
cd ../
:: build Module CalcStarter.exe  
cd ./CalcStarter 
"g++.exe" -static -O2 -Wall -std=c++14 main.cpp -o CalcStarter 2> !logg++.txt
copy CalcStarter.exe  %PATHTOCALCFOLDER%       
cd ../
