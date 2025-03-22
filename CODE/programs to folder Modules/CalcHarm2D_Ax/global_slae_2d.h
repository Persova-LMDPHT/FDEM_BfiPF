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
 *  This file contains headers for some basic routines for assembling the finite element matrix and vector of the right hand side
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin  
 *  Novosibirsk State Technical University,                                                                      
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                            
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                             
 *  Version 2.0 December 10, 2024                                                                                 
*/

#pragma once

class T_Global_SLAE_2d
{
public:
	long   *ig, *jg;
	double *di, *ggl, *ggu, *di_b, *ggl_b, *ggu_b;
	double *pr;
	long n_elem;
	long n_nodes;
	long n_nodes_c;
	long *bound_nodes;
	short *bound_mtrs;
	long n_bound_nodes;
	long n;
	long ig_n_1;
	long *node_is_bound;
	long (*nvtr)[4];
	double (*xy)[2];
	long *nvkat;
	long n_of_materials;
	double *sigma;
	double *dpr;
	double *mu;
	double omega;
	double *sigma0;
	long n_1d;
	double *coords_1d;
	double *sin_1d;
	double *cos_1d;
	bool type_of_task;
	long *type_of_rect;
	int alpha;
	bool is_RZ;
	double currentDeltaFunc;
	int npls;
	T_Global_SLAE_2d(long *ig, long *jg, long n_elem, long n_nodes, long n_bound_nodes,
		long (*nvtr)[4], double (*xy)[2], long *bound_nodes, short *bound_mtrs,
		long n_nodes_c,long n_materials, double *mu2d, long *nvkat,
		double nu, double *sigma2d, double *dpr2d, int alpha, long *type_of_rect,
		double *sigma0, long n_1d, double *coords_1d, double *sin_1d, double *cos_1d,
		double currentDeltaFunc, int _npls);
	~T_Global_SLAE_2d();
	// Assembling the global matrix of the finite element SLAE
	void Assembling_With_T_Mapping(bool Ax=false);
};
