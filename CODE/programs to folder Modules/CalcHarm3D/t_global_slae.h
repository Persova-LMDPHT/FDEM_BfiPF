/*
 * GENERAL REMARKS
 * 
 *  This code is freely available under the following conditions:
 *  
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *  
 *                          
 *  This file contains headers for some basic routines for assembling the finite element matrix and vector of the right hand side.
 * 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/t_global_slae.h 
 *  Novosibirsk State Technical University,                                                                  
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                        
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                         
 *
*/ 

#pragma once
#include "T_Brick.h"

class VecBlockSLAE
{
public:
	long   *ig, *jg;
	double *di_block;
	double *gg_block;
	long *idi;
	long *ijg;
	double *pr;
	long n_elem;
	long n_edges;
	long n_edges_c;
	long n;
	long nb;
	long ig_n_1;
	long (*nver)[14];
	long (*ed)[25];
	double (*xyz)[3];
	long (*edges)[2];
	long n_of_materials;
	long *nvkat;
	double *dpr3d;
	double *dpr0;
	double *sigma3d;
	double *sigma0;
	double *mu3d;
	double *mu0;
	double omega;
	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;
	int tasktype;
	VecBlockSLAE(long *ig, long *jg, long *idi, long *ijg,
		long n_elem, long n_edges, long n_edges_c,	
		double (*xyz)[3], long (*nver)[14], long (*ed)[25], long (*edges)[2], 
		long *nvkat, double nu, double *mu3d, double *mu0, double *sigma3d, double *sigma0, double *dpr3d, double *dpr0,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d,int npls);
	~VecBlockSLAE();
	// Assembling the global matrix and vector of the right side of the finite element SLAE
	void AsmBlockSLAE(Vec_Prep_Data *d);
	// Assembling the global matrix of the finite element SLAE
	void AsmBlockSLAE_MATRIX(Vec_Prep_Data *d);
	// Assembling the global vector of the right-hand side of the finite element SLAE
	void AsmBlockSLAE_PR(Vec_Prep_Data *d);
	void AsmBlockEnorm(Vec_Prep_Data *d,double *enorm);
	void Add_to_di_block(T_Brick *LocElem, int i, int j, long nBlock, double mult);
	void Add_to_gg_block(T_Brick *LocElem, int i, int j, long nBlock, double mult);
	void Add_to_pr_block(T_Brick *LocElem, int i, long nBlock, double mult);
	int npls;
};
