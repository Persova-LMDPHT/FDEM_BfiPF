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
 *  This file contains headers for working with edge mesh in 3D VFEM
 *  
 *  Written by D.Sc. Denis V. Vagin                                           
 *  Novosibirsk State Technical University,                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                         
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)          
 * Version 2.0 December 10, 2024                                              
 *  
*/

#pragma once

class T_Mapping_Vec
{
public: 
	T_Mapping_Vec();
	T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar);
	~T_Mapping_Vec();

	int n_c; 
	int n_dc;
	int n; 

	int kuzlov;
	int kpar; 
	int (*nver)[14];
	double (*xyz)[3];

	int (*edges)[2]; 
	int (*ed)[25]; 
};
