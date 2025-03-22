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
#include "For_Solvers.h"

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
void T_Global_SLAE_2d::Asm_B_C_Output(int *ig, int *jg, int n_elem,
											 int (*nvtr)[4], double (*xy)[2],
											 int n_nodes_c, double *mu2d, int *nvkat,
											 double *sigma2d, double *sigma2d_Z, double *dpr2d, double omega, double *w_re, double *w_im, int npls,int kuzlov)
{
	int i, j, k, m, ipls;
	int it, jt, i_mu, j_nu;
	int ii, jj;
	int i_glob, j_glob;
	int n = n_nodes_c;
	int ig_n_1 = ig[n];

	cout << "Asm_B_C_for_Output_Er...\n";
	logfile << "Asm_B_C_for_Output_Er...\n";

	if (di_b==NULL)
		if((di_b = new double[n*2])==NULL) Memory_allocation_error("di_b", "Asm_B_C_for_Output_Er...");
	if (ggl_b==NULL)
		if((ggl_b = new double[ig_n_1*2])==NULL) Memory_allocation_error("ggl_b", "Asm_B_C_for_Output_Er...");
	if (ggu_b==NULL)
		if((ggu_b = new double[ig_n_1*2])==NULL) Memory_allocation_error("ggu_b", "Asm_B_C_for_Output_Er...");
	if (pr==NULL)
		if((pr = new double[n*2*npls])==NULL) Memory_allocation_error("pr", "Asm_B_C_for_Output_Er...");

	for(i=0; i<n*2; i++)
	{
		di_b[i] = 0.0;
	}

	for(i=0; i<n*2*npls; i++)
	{
		pr[i] = 0.0;
	}

	for(i=0; i<ig_n_1*2; i++) 
	{
		ggl_b[i] = 0.0; 
		ggu_b[i] = 0.0; 
	}

	for(i=0; i<n_elem; i++)
	{
		double lmdr[4][4],lmdz[4][4],lmij,vdsig[4],d2wdrdzsm[2],vvdsig[4];

		double r0 = xy[nvtr[i][0]][0];
		double r1 = xy[nvtr[i][3]][0];
		double z0 = xy[nvtr[i][0]][1];
		double z1 = xy[nvtr[i][3]][1];
		double mu = mu2d[nvkat[i]];
		double sigma = sigma2d[nvkat[i]];
		double sigmaZ = sigma2d_Z[nvkat[i]];
		double dpr = dpr2d[nvkat[i]];

		LocalMatrixRZ L(r0, r1, z0, z1, mu, sigma, sigmaZ, 0);

		L.CalcLocalMatrixC();
		GetStiffnessMatrixDr(r0, z1, r1, z0, lmdr);
		GetStiffnessMatrixDz(r0, z1, r1, z0, lmdz);

		vvdsig[0]=0.0;
		vvdsig[1]=0.0;
		vvdsig[2]=0.0;
		vvdsig[3]=0.0;

		vdsig[0]=((sigmaZ-sigma)/(mu*sigma))*(-r1*r1-r1*r0+r0*r0*2.0)/6.0;
		vdsig[1]=((sigmaZ-sigma)/(mu*sigma))*(r1*r0-r1*r1*2.0+r0*r0)/6.0;
		vdsig[2]=((sigmaZ-sigma)/(mu*sigma))*(r1*r1+r1*r0-r0*r0*2.0)/6.0;
		vdsig[3]=((sigmaZ-sigma)/(mu*sigma))*(r1*r1*2.0-r1*r0-r0*r0)/6.0;

		for(j=0;j<4;j++)
		{
			lmij=lmdr[j][j]+lmdz[j][j]*(sigmaZ/sigma);

			ii = nvtr[i][j];
			if(ii < n_nodes_c)
			{
				i_glob = ii;
				di_b[i_glob*2]   += lmij/mu-L.c[j][j]*dpr*omega*omega;
				di_b[i_glob*2+1] += L.c[j][j]*sigmaZ*omega;
			}

			for(k=0;k<4;k++)
			{
				lmij=lmdr[j][k]+lmdz[j][k]*(sigmaZ/sigma);

				jj = nvtr[i][k];

				i_glob = nvtr[i][j];
				j_glob = nvtr[i][k];
				if(j_glob < i_glob)
				{
					for(m=ig[i_glob]; m<=ig[i_glob+1]-1; m++)
					{
						if(jg[m]==j_glob)
						{
							ggl_b[m*2]   += lmij/mu-L.c[j][k]*dpr*omega*omega;
							ggl_b[m*2+1] += L.c[j][k]*sigmaZ*omega;

							ggu_b[m*2]   += lmij/mu-L.c[k][j]*dpr*omega*omega;
							ggu_b[m*2+1] += L.c[k][j]*sigmaZ*omega;
						}
					}
				}
			}

			for(ipls=0;ipls<npls;ipls++)
			{
				d2wdrdzsm[0] = ((w_re[nvtr[i][3]+ipls*kuzlov] - w_re[nvtr[i][2]+ipls*kuzlov])/(r1-r0)+
					(w_re[nvtr[i][1]+ipls*kuzlov] - w_re[nvtr[i][0]+ipls*kuzlov])/(r1-r0))/2.0;
				d2wdrdzsm[1] = ((w_im[nvtr[i][3]+ipls*kuzlov] - w_im[nvtr[i][2]+ipls*kuzlov])/(r1-r0)+
					(w_im[nvtr[i][1]+ipls*kuzlov] - w_im[nvtr[i][0]+ipls*kuzlov])/(r1-r0))/2.0;

				d2wdrdzsm[0]*=vdsig[j];
				d2wdrdzsm[1]*=vdsig[j];

				i_glob = ii;
				pr[i_glob*2+ipls*n*2]   += d2wdrdzsm[0];
				pr[i_glob*2+ipls*n*2+1] += d2wdrdzsm[1];
			}
		}
	}	
}

