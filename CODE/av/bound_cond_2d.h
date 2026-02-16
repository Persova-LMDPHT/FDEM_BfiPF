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
 *  Version 2.0 December 10, 2024                                                                                   
*/

#pragma once

class Bound_Cond
{
public:
	double *di;
	double *ggl;
	double *ggu;
	double *pr;
	int n_nodes_c;
	int n_bound_nodes;
	int *bound_nodes;
	short *bound_mtrs;
	int *ig;
	int *jg;

	Bound_Cond(double *di, double *ggl, double *ggu, double *pr, int *bound_nodes,
		int n_nodes_c, int n_bound_nodes, int *ig, int *jg);
	~Bound_Cond();

	void Set_bound_cond();
	// Taking into account the first boundary conditions
	void Set_bound_cond_harm(int npls);
	void Set_bound_cond_harm(double (*xy)[2]);
	void Set_bound_cond_harm2(int n_2kr, int (*edges_2kr)[2]);
	void Set_bound_cond_with_symmetrization(double (*xy)[2]);

	void Set2kr(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
		double *mu2d, double current, int n_2kr, int (*edges_2kr)[2]);

	void Set2krHarm(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
		double *mu2d, double current, int n_2kr, int (*edges_2kr)[2],int ipls);

	void Set2kr_3x3(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
		double *mu2d, double current, int n_2kr, int (*edges_2kr)[2]);

	void Set1kr(double (*xy)[2]);

	void Set1kr_big(double (*xy)[2], double big);

	void SwitchOffEquations(bool a1, bool a2, bool a3);
	void SwitchOffEquationsHarm(bool a1, bool a2, bool a3);

	void SetDeltaFunction(int elem1, int elem2, double current, int (*nvtr)[4], double (*xy)[2], double omega);
	void SetDeltaFunction2(int node, double current, int (*nvtr)[4], double (*xy)[2], double omega);
	void SetDeltaFunctionElliptic(int node, double current, int (*nvtr)[4], double (*xy)[2], double omega);
	void SetF(int elem, double current, int (*nvtr)[4], double (*xy)[2]);
};
