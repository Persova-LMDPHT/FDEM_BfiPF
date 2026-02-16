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
 *  This file contains code of some basic routines for assembling the finite element matrix and vector of the right hand side
 *  
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,                                                                      
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                            
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                             
 *  Version 2.0 December 10, 2024                                                                                  
*/

#include "stdafx.h"
#include "global_slae_2d.h"
#include "in_out.h"
#include "For_Solvers.h"
#include "rect_local_matrix_2d.h"
#include "bound_cond_2d.h"

extern ofstream logfile;

extern int GetStiffnessMatrix(const double &r0, const double &z0,
					   const double &r3, const double &z3,
					   double m[4][4]);

extern int GetStiffnessMatrix0(const double &r0, const double &z0,
					   const double &r3, const double &z3,
					   double m[4][4]);

extern int GetMassMatrix(const double &r0, const double &z0,
				  const double &r3, const double &z3,
				  double m[4][4]);

T_Global_SLAE_2d::T_Global_SLAE_2d(long *ig, long *jg, long n_elem, long n_nodes, long n_bound_nodes,
								   long (*nvtr)[4], double (*xy)[2], long *bound_nodes, short *bound_mtrs,
								   long n_nodes_c,long n_materials, double *mu2d, long *nvkat,
								   double nu, double *sigma2d, double *dpr2d, int alpha, long *type_of_rect,
								   double *sigma0, long n_1d,
								   double *coords_1d, double *sin_1d, double *cos_1d,
								   double currentDeltaFunc, int _npls)
{
	this->ig = ig;
	this->jg = jg;
	this->n_elem = n_elem;

	this->n_nodes = n_nodes;
	this->n_bound_nodes = n_bound_nodes;

	this->nvtr = nvtr;
	this->xy = xy;
	this->bound_nodes = bound_nodes;
	this->bound_mtrs = bound_mtrs;

	this->n_nodes_c = n_nodes_c;
	this->n = n_nodes_c;
	this->ig_n_1 = ig[n];

	this->n_of_materials = n_materials;
	this->mu = mu2d;
	this->nvkat = nvkat;
	this->omega = nu*2.0*3.1415926535;
	this->sigma = sigma2d;
	this->dpr = dpr2d;

	this->alpha = alpha;

	this->type_of_rect = type_of_rect;

	this->sigma0 = sigma0;
	this->n_1d = n_1d;
	this->coords_1d = coords_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;

	pr=NULL;
	di=NULL;
	ggl=NULL;
	ggu=NULL;

	di_b = NULL;
	ggl_b = NULL;
	ggu_b = NULL;

	is_RZ = true;
	this->currentDeltaFunc = currentDeltaFunc;

	npls=_npls;
}

T_Global_SLAE_2d::~T_Global_SLAE_2d()
{
	if(pr) {delete [] pr; pr=NULL;}
	if(di) {delete [] di; di=NULL;}
	if(ggl) {delete [] ggl; ggl=NULL;}
	if(ggu) {delete [] ggu; ggu=NULL;}
	if (di_b) { delete[] di_b; di_b = NULL; }
	if (ggl_b) { delete[] ggl_b; ggl_b = NULL; }
	if (ggu_b) { delete[] ggu_b; ggu_b = NULL; }
}
//------------------------------------------------------ 
// Assembling the global matrix of the finite element SLAE
//------------------------------------------------------ 
void T_Global_SLAE_2d::Assembling_With_T_Mapping(bool Ax)
{
	long i, j, k, m;
	long k2, j2;
	long it, jt, i_mu, j_nu;
	long ii, jj;
	long i_glob, j_glob;
	long loc_str, loc_col;
	const double dpr0=8.84194128288307421e-12; 

	cout << "Assembling global SLAE for MTZ using T_mapping...\n";
	logfile << "Assembling global SLAE for MTZ using T_mapping...\n";

	if (di_b == NULL)
		if ((di_b = new double[n*2]) == NULL) Memory_allocation_error("di_b", "Assembling_With_T_Mapping");
	if (ggl_b == NULL)
		if ((ggl_b = new double[ig_n_1*2]) == NULL) Memory_allocation_error("ggl_b", "Assembling_With_T_Mapping");
	if (ggu_b == NULL)
		if ((ggu_b = new double[ig_n_1*2]) == NULL) Memory_allocation_error("ggu_b", "Assembling_With_T_Mapping");
	if (pr == NULL)
		if ((pr = new double[n*2*npls]) == NULL) Memory_allocation_error("pr", "Assembling_With_T_Mapping");

	for (i=0;i<n*2;i++){di_b[i]=0.0;}
	for (i=0;i<n*2*npls;i++){pr[i]=0.0;}
	for (i=0;i<ig_n_1*2;i++){ggl_b[i]=ggu_b[i]=0.0;}

	for(i=0; i<n_elem; i++)
	{
		Rect_Local_Matrix L;

		if (is_RZ)
		{
			L = Rect_Local_Matrix(i, xy, nvtr);

			(Ax) ?
				GetStiffnessMatrix0(L.x[0], L.y[0], L.x[3], L.y[3], L.b)
				:
			GetStiffnessMatrix(L.x[0], L.y[0], L.x[3], L.y[3], L.b);
		
			GetMassMatrix(L.x[0], L.y[0], L.x[3], L.y[3], L.c);

			double _mu = mu[nvkat[i]];
			double _sigma = sigma[nvkat[i]];
			double _dpr = dpr[nvkat[i]]*dpr0;
			
			double lambda_1,lambda_2;
			lambda_1=1.0/_mu;
			lambda_2=0.0;

			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					L.a[i*2][j*2] = L.b[i][j]*lambda_1-L.c[i][j]*_dpr*omega*omega;
					L.a[i*2][j*2+1] = -L.c[i][j]*_sigma*omega+L.b[i][j]*lambda_2;
					L.a[i*2+1][j*2] = L.c[i][j]*_sigma*omega-L.b[i][j]*lambda_2;
					L.a[i*2+1][j*2+1] = L.b[i][j]*lambda_1-L.c[i][j]*_dpr*omega*omega;
				}
			}
		} 
		else
		{
			L = Rect_Local_Matrix(i, xy, nvtr, type_of_rect, mu,
				sigma, omega, nvkat, alpha, sigma0, n_1d, coords_1d, sin_1d, cos_1d);
			L.Calc_local_matrix();
		}

		for(j=0;j<4;j++)
		{
			loc_str=j*2;

			i_glob = nvtr[i][j];
			di_b[i_glob*2] += L.a[loc_str][loc_str];
			di_b[i_glob*2+1] += L.a[loc_str+1][loc_str];
	
			for(k=0;k<4;k++)
			{
				loc_col=k*2;

				j_glob = nvtr[i][k];
				if (j_glob < i_glob)
				{
					for (m = ig[i_glob]; m <= ig[i_glob + 1] - 1; m++)
					{
						if (jg[m] == j_glob)
						{
							ggl_b[m*2] += L.a[loc_str][loc_col];
							ggl_b[m*2+1] += L.a[loc_str+1][loc_col];

							ggu_b[m*2] += L.a[loc_col][loc_str];
							ggu_b[m*2+1] += L.a[loc_col+1][loc_str];
						}
					}
				}
			}
		}
	}
}
