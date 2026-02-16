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
	void Set_bound_cond_harm();
	void Set_bound_cond_harm(double (*xy)[2]);
	void Set_bound_cond_harm2(int n_2kr, int (*edges_2kr)[2]);
	void Set_bound_cond_with_symmetrization(double (*xy)[2]);

	// Taking into account the first boundary conditions
	void Set_bound_cond_harm_output(double (*u_bnd)[2],int npls);

	void SwitchOffEquations(bool a1, bool a2, bool a3);
	void SwitchOffEquationsHarm(bool a1, bool a2, bool a3);

	void CalcLocalMatrixUp(double (*m)[4], double hr, double hz, double rk, double gamma);
	void CalcLocalMatrixDown(double (*m)[4], double hr, double hz, double rk, double gamma);
	void CalcLocalVectorUp(double *b, double hr, double rk, double gamma_re, double gamma_im);
	void CalcLocalVectorErUp(double *b, double hr, double rk, double er_re, double er_im);
	void CalcLocalVectorErDown(double *b, double hr, double rk, double er_re, double er_im);
	void CalcLocalVectorDown(double *b, double hr, double rk, double gamma_re, double gamma_im);
// Addition to the matrix and the right-hand side of the term with a jump in sigma at the separation of the media
	int SetSigma(int nelmJ, int *elmJ, double (*sigJ)[2], double (*epsJ)[2], int *lowUpElm, double k, int (*nvtr)[4],
		double (*xy)[2], double *w_re, double *w_im, int npls, int kuzlov, int *nvkat, double *sigma2d, double *sigma2d_Z, double *dpr2d, double omega);
	void AddLocalMatrix(double (*m)[4], int (*nvtr)[4], int elm, double gamma_re,double gamma_im);
	void AddLocalVector(double *v, int (*nvtr)[4], int elm, int ipls);
};
