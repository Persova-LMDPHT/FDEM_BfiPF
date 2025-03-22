/*
 * GENERAL REMARKS
 *  
 *  This code is freely available under the following conditions:
 *  
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *  
 *                       FDEM_BfiPF   			
 *  This file contains headers for reading and storing edge mesh in 3D VFEM
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                            
 *  Novosibirsk State Technical University,                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                          
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)           
 *  Version 2.0 December 10, 2024                                              
*/

#pragma once
#include "T_Mapping.h"
#include "Tensor.h"

class Vec_Prep_Data
{
public:
	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data();
	int Read_mesh_for_nonstat_problem(char *pointres_fname);

	int AnomalType,fda;

	int Read_3dmeshregular(int interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	int maxiter;
	double eps;
	int n_materials;
	double *mu3d;
	double *mu0;
	int n_pointresB;
	int n_pointresE;
	double (*pointresB)[3];
	double (*pointresE)[3];
	double *sigma3d;
	double *sigma0;
	Tensor *sigmaTensor;
	Tensor *sigma0Tensor;
	vector<int> isTensor;
	double *dpr3d;
	double *dpr0;
	int kuzlov;
	int kpar;
	int kt1;
	int *l13d;
	int (*nver)[14];
	int *nvkat;
	double (*xyz)[3];
	int n_layers_1d;
	double *layers_1d;
	double *sigma_1d;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	int n_mesh_regular_x;
	int n_mesh_regular_y;
	double *mesh_regular_x;
	double *mesh_regular_y;

	int ntime;
	double *time;

 	int norvect;
 	double nu;
 	int alfa;
 	double *usin;
 	double *ucos;
 	double *z_1d;
 	int n_1d;

	T_Mapping_Vec *tmap;
	int tasktype;

	int npr, nfreq;

	int nobj;
	vector<double> Xobj[2],Yobj[2],Zobj[2];
	vector<int> Mobj;
};
