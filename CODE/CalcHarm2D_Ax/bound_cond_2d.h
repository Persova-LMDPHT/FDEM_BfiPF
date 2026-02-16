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
 *  The file contains headers for the functions of the class "Bound_cond"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin     
 *  Novosibirsk State Technical University,                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                        
 * Version 2.0 December 10, 2024                                                              
*/

#pragma once

class Bound_Cond
{
public:
	double *di;
	double *ggl;
	double *ggu;
	double *pr;
	long n_nodes_c;
	long n_bound_nodes;
	long *bound_nodes;
	short *bound_mtrs;
	long *ig;
	long *jg;
	Bound_Cond(double *di, double *ggl, double *ggu, double *pr, long *bound_nodes, short *bound_mtrs,
		long n_nodes_c, long n_bound_nodes, long *ig, long *jg);
	~Bound_Cond();
	void Set_bound_cond_for_Ay_problem(double (*xy)[2], long (*nvtr)[4], long n_elem);
	// Taking into account the first boundary conditions
	void Set_bound_cond_for_Hy_problem_with_symmetrization(int npls, double (*xy)[2], double cd, double *di_b, double *ggl_b, double *ggu_b);
};
