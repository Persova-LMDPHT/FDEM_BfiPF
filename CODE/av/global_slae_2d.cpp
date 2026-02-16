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
#include "Local_matrix_rz.h"

extern ofstream logfile;

T_Global_SLAE_2d::T_Global_SLAE_2d()
{
	pr=NULL;
	di_b=NULL;
	ggl_b=NULL;
	di_c=NULL;
	ggl_c=NULL;
	ggu_b=NULL;
}

void T_Global_SLAE_2d::Clear()
{
	if(pr) {delete [] pr; pr=NULL;}
	if(di_b) {delete [] di_b; di_b=NULL;}
	if(ggl_b) {delete [] ggl_b; ggl_b=NULL;}
	if(ggu_b) {delete [] ggu_b; ggu_b=NULL;}
	if(di_c) {delete [] di_c; di_c=NULL;}
	if(ggl_c) {delete [] ggl_c; ggl_c=NULL;}
}

T_Global_SLAE_2d::~T_Global_SLAE_2d()
{
	Clear();
}
//------------------------------------------------------ 
// Assembling the global matrix of the finite element SLAE
//------------------------------------------------------ 
void T_Global_SLAE_2d::AsmHarmDivMatrix(int *ig, int *jg, int n_elem,
								  int (*nvtr)[4], double (*xy)[2], int n_nodes_c, double *mu2d, int *nvkat,
								  vector<double> &sigma_xy, vector<double> &sigma_z, double omega, double *pr_val)
{
	int i, j, k, m;
	int k2, j2;
	int it, jt, i_mu, j_nu;
	int ii, jj, kk;
	int i_glob, j_glob;
	int loc_str, loc_col;
	int n = n_nodes_c*3;
	int ig_n_1 = ig[n];
	double pr_val_r_re[4];
	double pr_val_r_im[4];

	cout << "Assembling global SLAE using T_mapping...\n";
	logfile << "Assembling global SLAE using T_mapping...\n";

	if (di_b==NULL)
		if((di_b = new double[n*2])==NULL) Memory_allocation_error("di_b", "T_Global_SLAE_2d::Asm_B_C_Pr");
	if (ggl_b==NULL)
		if((ggl_b = new double[ig_n_1*2])==NULL) Memory_allocation_error("ggl_b", "T_Global_SLAE_2d::Asm_B_C_Pr");
	if (ggu_b==NULL)
		if((ggu_b = new double[ig_n_1*2])==NULL) Memory_allocation_error("ggu_b", "T_Global_SLAE_2d::Asm_B_C_Pr");

	for(i=0; i<n*2; i++)
	{
		di_b[i] = 0.0;
	}

	for(i=0; i<ig_n_1*2; i++) 
	{
		ggl_b[i] = 0.0; 
		ggu_b[i] = 0.0; 
	}

	for(i=0; i<n_elem; i++)
	{
		double r0 = xy[nvtr[i][0]][0];
		double r1 = xy[nvtr[i][3]][0];
		double z0 = xy[nvtr[i][0]][1];
		double z1 = xy[nvtr[i][3]][1];
		double mu = mu2d[nvkat[i]];
		double sig_xy = sigma_xy[nvkat[i]];
		double sig_z = sigma_z[nvkat[i]];

		LocalMatrixRZ L(r0, r1, z0, z1, mu, 0, omega);
		L.CalcLocalMatrixHarm(omega, sig_xy, sig_z);

		for (j=0; j<4; j++)
		{
			pr_val_r_re[j] = pr_val[nvtr[i][j]*2];
			pr_val_r_im[j] = pr_val[nvtr[i][j]*2+1];
		}
		L.CalcLocalVectorHarmDiv(pr_val_r_re, pr_val_r_im, sig_xy, sig_z);

		for(j=0; j<4; j++)
		{
			ii = nvtr[i][j];
			for(j2=0; j2<3; j2++)
			{
				loc_str = j*3 + j2;

				i_glob = nvtr[i][j]*3 + j2;
				for (kk=0; kk<2; kk++)
				{
					di_b[i_glob*2+kk] += L.ah[loc_str][loc_str][kk];
				}

				for(k=0; k<4; k++)
				{
					jj = nvtr[i][k];
					for(k2=0; k2<3; k2++)
					{
						loc_col = k*3 + k2;
						i_glob = nvtr[i][j]*3 + j2;
						j_glob = nvtr[i][k]*3 + k2;
						if (j_glob < i_glob)
						{
							for(m=ig[i_glob]; m<=ig[i_glob+1]-1; m++)
							{
								if(jg[m]==j_glob)
								{
									for (kk=0; kk<2; kk++)
									{
										ggl_b[m*2+kk] += L.ah[loc_str][loc_col][kk];
										ggu_b[m*2+kk] += L.ah[loc_col][loc_str][kk];
									}
									break;
								}
							}
						}
					}
				}
			}
		}
	}
}
//------------------------------------------------------ 
// Assembling the global vector of the right-hand side of the finite element SLAE
//------------------------------------------------------ 
void T_Global_SLAE_2d::AsmHarmDivPr(int *ig, int *jg, int n_elem,
							   int (*nvtr)[4], double (*xy)[2],int n_nodes_c, int kuzlov,  double *mu2d, int *nvkat,
							   vector<double> &sigma_xy, vector<double> &sigma_z, double omega, double *pr_val,int npls)
{
	int i, j, k, m;
	int k2, j2;
	int it, jt, i_mu, j_nu;
	int ii, jj;
	int i_glob, j_glob;
	int loc_str, loc_col;
	int n = n_nodes_c*3;
	int ig_n_1 = ig[n];
	double pr_val_r_re[4];
	double pr_val_r_im[4];
	int ipls;

	cout << "Assembling global SLAE using T_mapping...\n";
	logfile << "Assembling global SLAE using T_mapping...\n";

	if (pr==NULL)
		if((pr = new double[n*2*npls])==NULL) Memory_allocation_error("pr", "T_Global_SLAE_2d::Asm_B_C_Pr");

	for(i=0; i<n*2*npls; i++)
	{
		pr[i] = 0.0;
	}

	for(i=0; i<n_elem; i++)
	{
		double r0 = xy[nvtr[i][0]][0];
		double r1 = xy[nvtr[i][3]][0];
		double z0 = xy[nvtr[i][0]][1];
		double z1 = xy[nvtr[i][3]][1];
		double mu = mu2d[nvkat[i]];
		double sig_xy = sigma_xy[nvkat[i]];
		double sig_z = sigma_z[nvkat[i]];

		LocalMatrixRZ L(r0, r1, z0, z1, mu, 0, omega);
		L.CalcLocalMatrixHarm(omega, sig_xy, sig_z);

		for(ipls=0;ipls<npls;ipls++)
		{
			for (j=0; j<4; j++)
			{
				pr_val_r_re[j] = pr_val[nvtr[i][j]*2+ipls*2*kuzlov];
				pr_val_r_im[j] = pr_val[nvtr[i][j]*2+1+ipls*2*kuzlov];
			}
			L.CalcLocalVectorHarmDiv(pr_val_r_re, pr_val_r_im, sig_xy, sig_z);

			for(j=0; j<4; j++)
			{
				ii = nvtr[i][j];
				for(j2=0; j2<3; j2++)
				{
					loc_str = j*3 + j2;
					i_glob = nvtr[i][j]*3 + j2;

					for (k=0; k<2; k++)
					{
						pr[i_glob*2+k+ipls*2*n]   += L.gh[loc_str][k];
					}
				}
			}
		}
	}
}
