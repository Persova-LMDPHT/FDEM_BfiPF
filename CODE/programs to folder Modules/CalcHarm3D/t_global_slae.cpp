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
 *  This file contains code of some basic routines for assembling the finite element matrix and vector of the right hand side.
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/t_global_slae.cpp
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *
*/

#include "stdafx.h"
#include "t_global_slae.h"

extern ofstream logfile;

VecBlockSLAE::VecBlockSLAE(long *ig, long *jg, long *idi, long *ijg, 
			  long n_elem, long n_edges, long n_edges_c,	
			  double (*xyz)[3], long (*nver)[14], long (*ed)[25], long (*edges)[2], 
			  long *nvkat, double nu, double *mu3d, double *mu0, double *sigma3d, double *sigma0, double *dpr3d, double *dpr0,
			  long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d,int _npls)
{
	this->ig = ig;
	this->jg = jg;
	this->idi = idi;
	this->ijg = ijg;

	this->n_elem = n_elem;
	this->n_edges = n_edges; 
	this->n_edges_c = n_edges_c; 

	this->xyz = xyz;
	this->nver = nver;
	this->ed = ed;
	this->edges = edges;

	this->nvkat = nvkat;
	this->omega = 2.0*PI*nu;
	this->mu3d = mu3d;
	this->mu0 = mu0;
	this->sigma3d = sigma3d;
	this->sigma0 = sigma0;
	this->dpr3d = dpr3d;
	this->dpr0 = dpr0;

	this->alpha = alpha;
	this->n_1d = n_1d;
	this->z_1d = z_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;

	this->nb = n_edges_c;
	this->n = n_edges_c*2;

	this->ig_n_1 = this->ig[this->n_edges_c];

	npls=_npls;

	di_block = NULL;
	gg_block = NULL;
	pr = NULL;
	if((pr = new double[this->n*npls])==0) Memory_allocation_error("pr", "VecBlockSLAE::VecBlockSLAE");
	if((di_block = new double[idi[n_edges_c]])==0) Memory_allocation_error("di_block", "VecBlockSLAE::VecBlockSLAE");
	if((gg_block = new double[ijg[ig_n_1]])==0) Memory_allocation_error("gg_block", "VecBlockSLAE::VecBlockSLAE");

	tasktype=0;
}

VecBlockSLAE::~VecBlockSLAE()
{
	if(pr)  { delete [] pr; pr=NULL; }
	if(gg_block) {delete [] gg_block; gg_block=NULL;}
	if(di_block) {delete [] di_block; gg_block=NULL;}
}
//-----------------------------------------------------------
// Assembling the global matrix and vector of the right side of the finite element SLAE
//-----------------------------------------------------------
void VecBlockSLAE::AsmBlockSLAE(Vec_Prep_Data *d)
{
	long i, j, k, m, it, jt, i_mu, j_nu, ipls, st, l;
	long ii, jj;
	for(i=0; i<=idi[nb]-1; i++)
		di_block[i] = 0;

	for(i=0; i<=ijg[ig_n_1]-1; i++)
		gg_block[i] = 0;

	for(i=0; i<n*npls; i++)
		pr[i] = 0.0;

	for(i=0; i<n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
			alpha, n_1d, z_1d, sin_1d, cos_1d);

		L.d=d;
		L.tasktype=tasktype;

		L.sigmaTensor=d->sigma3dTensor[nvkat[i]];
		L.sigmaTensor0=d->sigma2dTensor[nvkat[i]];
		L.dprTensor=dpr3d[nvkat[i]];
		L.dprTensor0=dpr0[nvkat[i]];

		L.GetRotationMatrix();
		L.RotateTensors();

		L.Calc_ss0();

		L.mu0=mu0[nvkat[i]];
		L.mu=mu3d[nvkat[i]];			
		L.sigma0=sigma0[nvkat[i]];
		L.sigma=sigma3d[nvkat[i]];
		L.dpr0=dpr0[nvkat[i]];
		L.dpr=dpr3d[nvkat[i]];

		L.ipls=0;
		L.Calc_block_local_matrix_and_vector(true,false);

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii >= n_edges_c)
			{
			}
			else
			{  
				Add_to_di_block(&L, j,j, ii, 1.0); 
			}

			for(k=0; k<12; k++)
			{
				jj = ed[i][k];

				if(jj < ii)
				{
					for(m=ig[ii]; m<=ig[ii+1]-1; m++)
					{
						if(jg[m]==jj)
						{
							Add_to_gg_block(&L, j,k, m, 1.0);
							break;
						}
					}
				}
			}
		}

		if(!(d->fdirect))
		{
			L.n_edges=n_edges;
			for(ipls=0;ipls<npls;ipls++)
			{
				L.ipls=ipls;
				L.Calc_block_local_matrix_and_vector(false,true);
				for(j=0; j<12; j++)
				{
					ii = ed[i][j];
					Add_to_pr_block(&L, j, ii+ipls*n_edges_c, 1.0);
				}
			}
		}
	}

	if(d->fdirect)
	{
		double val[3],bf[3],in[3];
		for(ipls=0;ipls<npls;ipls++)
		{
			st=(int)d->srs[ipls].size();
			for(l=0;l<st;l++)
			{
				loc_source &ls=d->srs[ipls][l];
				i=ls.elem;

				T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
					alpha, n_1d, z_1d, sin_1d, cos_1d);

				L.n_edges=n_edges;
				L.ipls=ipls;

				val[0]=ls.d_ksi*ls.len;
				val[1]=ls.d_eta*ls.len;
				val[2]=ls.d_dzeta*ls.len;

				in[0] = 2.0 * ls.ksi - 1.0;
				in[1]=2.0*ls.eta-1.0;
				in[2]=2.0*ls.dzeta-1.0;

				for(j=0; j<12; j++)
				{
					int nn[2];
					double ss[3];

					nn[0]=nver[i][REG_EDGES[j][0]];
					nn[1]=nver[i][REG_EDGES[j][1]];

					ss[0]=xyz[nn[1]][0]-xyz[nn[0]][0];
					ss[1]=xyz[nn[1]][1]-xyz[nn[0]][1];
					ss[2]=xyz[nn[1]][2]-xyz[nn[0]][2];

					L.Basis_func_on_reference_vec_par(j,in,bf);
					L.g_harm[j*2]=Scal(val,bf,3)*2.0/Norm_Euclid(ss,3);
					L.g_harm[j*2+1]=0.0;

					ii = ed[i][j];
					Add_to_pr_block(&L, j, ii+ipls*n_edges_c, 1.0);
				}
			}
		}
	}
}
//-----------------------------------------------------------
// Assembling the global matrix of the finite element SLAE 
//----------------------------------------------------------- 
void VecBlockSLAE::AsmBlockSLAE_MATRIX(Vec_Prep_Data *d)
{
	long i, j, k, m, it, jt, i_mu, j_nu, ipls, st, l;
	long ii, jj;

	for(i=0; i<=idi[nb]-1; i++)
		di_block[i] = 0;

	for(i=0; i<=ijg[ig_n_1]-1; i++)
		gg_block[i] = 0;

	for(i=0; i<n*npls; i++)
		pr[i] = 0.0;

	for(i=0; i<n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
			alpha, n_1d, z_1d, sin_1d, cos_1d);

		L.d=d;
		L.tasktype=tasktype;

		L.sigmaTensor=d->sigma3dTensor[nvkat[i]];
		L.sigmaTensor0=d->sigma2dTensor[nvkat[i]];
		L.dprTensor=dpr3d[nvkat[i]];
		L.dprTensor0=dpr0[nvkat[i]];

		L.GetRotationMatrix();
		L.RotateTensors();

		L.Calc_ss0();

		L.mu0=mu0[nvkat[i]];
		L.mu=mu3d[nvkat[i]];			
		L.sigma0=sigma0[nvkat[i]];
		L.sigma=sigma3d[nvkat[i]];
		L.dpr0=dpr0[nvkat[i]];
		L.dpr=dpr3d[nvkat[i]];

		L.ipls=0;
		L.Calc_block_local_matrix_and_vector(true,false);

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii >= n_edges_c)
			{
			}
			else
			{  
				Add_to_di_block(&L, j,j, ii, 1.0); 
			}


			for(k=0; k<12; k++)
			{
				jj = ed[i][k];

				if(jj < ii)
				{
					for(m=ig[ii]; m<=ig[ii+1]-1; m++)
					{
						if(jg[m]==jj)
						{
							Add_to_gg_block(&L, j,k, m, 1.0);
							break;
						}
					}
				}
			}
		}
	}
}
//-----------------------------------------------------------
// Assembling the global vector of the right-hand side of the finite element SLAE
//-----------------------------------------------------------
void VecBlockSLAE::AsmBlockSLAE_PR(Vec_Prep_Data *d)
{
	long i, j, k, m, it, jt, i_mu, j_nu, ipls, st, l;
	long ii, jj;

	for(i=0; i<n*npls; i++)
		pr[i] = 0.0;

	for(i=0; i<n_elem; i++)
	{
		T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
			alpha, n_1d, z_1d, sin_1d, cos_1d);

		L.d=d;
		L.tasktype=tasktype;

		L.sigmaTensor=d->sigma3dTensor[nvkat[i]];
		L.sigmaTensor0=d->sigma2dTensor[nvkat[i]];
		L.dprTensor=dpr3d[nvkat[i]];
		L.dprTensor0=dpr0[nvkat[i]];

		L.GetRotationMatrix();
		L.RotateTensors();

		L.Calc_ss0();

		L.mu0=mu0[nvkat[i]];
		L.mu=mu3d[nvkat[i]];			
		L.sigma0=sigma0[nvkat[i]];
		L.sigma=sigma3d[nvkat[i]];
		L.dpr0=dpr0[nvkat[i]];
		L.dpr=dpr3d[nvkat[i]];

		L.ipls=0;
		L.Calc_block_local_matrix_and_vector(true,false);

		if(!(d->fdirect))
		{
			L.n_edges=n_edges;
			for(ipls=0;ipls<npls;ipls++)
			{
				L.ipls=ipls;
				L.Calc_block_local_matrix_and_vector(false,true);
				for(j=0; j<12; j++)
				{
					ii = ed[i][j];
					Add_to_pr_block(&L, j, ii+ipls*n_edges_c, 1.0);
				}
			}
		}
	}

	if(d->fdirect)
	{
		double val[3],bf[3],in[3];
		for(ipls=0;ipls<npls;ipls++)
		{
			st=(int)d->srs[ipls].size();
			for(l=0;l<st;l++)
			{
				loc_source &ls=d->srs[ipls][l];
				i=ls.elem;

				T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
					alpha, n_1d, z_1d, sin_1d, cos_1d);

				L.n_edges=n_edges;
				L.ipls=ipls;

				val[0]=ls.d_ksi*ls.len;
				val[1]=ls.d_eta*ls.len;
				val[2]=ls.d_dzeta*ls.len;

				in[0]=2.0*ls.ksi-1.0;
				in[1]=2.0*ls.eta-1.0;
				in[2]=2.0*ls.dzeta-1.0;

				for(j=0; j<12; j++)
				{
					int nn[2];
					double ss[3];

					nn[0]=nver[i][REG_EDGES[j][0]];
					nn[1]=nver[i][REG_EDGES[j][1]];

					ss[0]=xyz[nn[1]][0]-xyz[nn[0]][0];
					ss[1]=xyz[nn[1]][1]-xyz[nn[0]][1];
					ss[2]=xyz[nn[1]][2]-xyz[nn[0]][2];

					L.Basis_func_on_reference_vec_par(j,in,bf);
					L.g_harm[j*2]=Scal(val,bf,3)*2.0/Norm_Euclid(ss,3);
					L.g_harm[j*2+1]=0.0;

					ii = ed[i][j];
					Add_to_pr_block(&L, j, ii+ipls*n_edges_c, 1.0);
				}
			}
		}
	}
}
void VecBlockSLAE::AsmBlockEnorm(Vec_Prep_Data *d,double *enorm)
{
	long i, j, k, m, it, jt, i_mu, j_nu, ipls;
	long ii, jj;

	if(!d->fdirect)
	{
		for(i=0; i<n_elem; i++)
		{
			T_Brick L(i, nver, ed, edges, xyz, nvkat, sigma3d, sigma0, mu3d, omega,
				alpha, n_1d, z_1d, sin_1d, cos_1d);

			L.d=d;
			L.tasktype=tasktype;
			L.n_edges=n_edges;

			for(ipls=0;ipls<npls;ipls++)
			{
				L.ipls=ipls;
				L.Calc_asin_acos_at_middle_of_edges();

				for(j=0; j<12; j++)
				{
					ii = ed[i][j]+n_edges*ipls;
					if(tasktype!=2)
					{
						enorm[ii*2]=L.acos0[j]*omega;
						enorm[ii*2+1]=-L.asin0[j]*omega;
					}
					else
					{
						enorm[ii*2]=L.asin0[j];
						enorm[ii*2+1]=L.acos0[j];
					}
				}
			}
		}
	}
}

void VecBlockSLAE::Add_to_di_block(T_Brick *LocElem, int i, int j, long nBlock, double mult)
{
	long beg = idi[nBlock];
	long size = idi[nBlock+1] - beg;

	if(size==1)
	{
		di_block[beg] += LocElem->b[i][j]*mult;  	
	}
	else
	{
		di_block[beg] += LocElem->b[i][j]*mult;
		di_block[beg+1] += LocElem->c[i][j]*mult;  	
	}
}

void VecBlockSLAE::Add_to_gg_block(T_Brick *LocElem, int i, int j, long nBlock, double mult)
{
	long beg = ijg[nBlock];
	long size = ijg[nBlock+1] - beg;

	if(size==1)
	{
		gg_block[beg] += LocElem->b[i][j]*mult;  	
	}
	else
	{
		gg_block[beg] += LocElem->b[i][j]*mult;
		gg_block[beg+1] += LocElem->c[i][j]*mult;  	
	}
}

void VecBlockSLAE::Add_to_pr_block(T_Brick *LocElem, int i, long nBlock, double mult)
{
	pr[nBlock*2]   += LocElem->g_harm[i*2]*mult;
	pr[nBlock*2+1] += LocElem->g_harm[i*2+1]*mult;
}
