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
	double *di_b, *ggl_b;
	double *di_c, *ggl_c;
	double *pr;
	double *ggu_b;

	T_Global_SLAE_2d();

	// Assembling the global matrix of the finite element SLAE
	void Asm_B_C_Output(int *ig, int *jg, int n_elem,
		int (*nvtr)[4], double (*xy)[2],
		int n_nodes_c,double *mu2d, int *nvkat,
		double *sigma2d, double *sigma2d_Z, double *dpr2d, double omega, double *w_re, double *w_im, int npls, int kuzlov);

	void Clear();

	~T_Global_SLAE_2d();
};
