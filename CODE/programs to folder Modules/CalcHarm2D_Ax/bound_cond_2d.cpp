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
 *  The file contains code of functions of the class "Bound_cond"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin     
 *  Novosibirsk State Technical University,                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                           
 * Version 2.0 December 10, 2024                                                             
*/

#include "stdafx.h"
#include "bound_cond_2d.h"
#include "For_Solvers.h"

Bound_Cond::Bound_Cond(double *di, double *ggl, double *ggu, double *pr,
					   long *bound_nodes, short *bound_mtrs, long n_nodes_c, long n_bound_nodes,
					   long *ig, long *jg)
{
	this->di = di;
	this->ggl = ggl;
	this->ggu = ggu;
	this->pr = pr;
	this->bound_nodes = bound_nodes;
	this->bound_mtrs = bound_mtrs;
	this->n_nodes_c = n_nodes_c;
	this->n_bound_nodes = n_bound_nodes;
	this->ig = ig;
	this->jg = jg;
}

Bound_Cond::~Bound_Cond()
{
}

void Bound_Cond::Set_bound_cond_for_Ay_problem(double (*xy)[2], long (*nvtr)[4], long n_elem)
{
	long i, j, node;
	bool *is_node_bound=NULL;

	is_node_bound = new bool[n_nodes_c];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "Bound_Cond::Set_bound_cond_for_Ay_problem");

	for(i=0; i<n_nodes_c; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = true;

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];

		for(j=0; j<2; j++)
		{
			long str = node*2+j;
			di[str] = 1.0;
			pr[str] = 0.0;

			for(long i1=ig[str]; i1<=ig[str+1]-1; i1++)
			{
				ggl[i1] = 0.0;
			}
		}
	}

	for(long i1=0; i1<n_nodes_c*2; i1++)
	{
		for(long j1=ig[i1]; j1<=ig[i1+1]-1; j1++)
		{
			if(is_node_bound[jg[j1]/2]==true)
			{
				ggu[j1] = 0.0;
			}
		}
	}
	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
//------------------------------------------------------ 
// Taking into account the first boundary conditions
//------------------------------------------------------ 
void Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization(int npls, double (*xy)[2], double cd,double *di_b, double *ggl_b, double *ggu_b)
{
	long i, j, k;
	long node, str;
	double u_re, u_im;
	double u_g;
	double *dirichlet=NULL;
	bool *is_dof_bound=NULL;
	int ipls;

	is_dof_bound = new bool[n_nodes_c];
	if(is_dof_bound == 0)
		Memory_allocation_error("is_dof_bound", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c; i++)
		is_dof_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
	{
		if (bound_mtrs&&(bound_mtrs[i]==2||bound_mtrs[i]==3))
			continue;
		is_dof_bound[bound_nodes[i]] = true;
	}

	dirichlet = new double[n_nodes_c*2];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_bound_nodes; i++)
	{
		node = bound_nodes[i];
		if (bound_mtrs&&(bound_mtrs[i]==2||bound_mtrs[i]==3))
			continue;
		if (bound_mtrs&&bound_mtrs[i]==1)
		{
			double rcoord=xy[node][0];
			dirichlet[node*2] = cd/(2*PI*rcoord);
		}
		else
			dirichlet[node*2] = 0.0;
		dirichlet[node*2+1] = 0.0; 
	}

	for (i = 0; i < n_bound_nodes; i++)
	{
		node = bound_nodes[i];

		if (is_dof_bound[node])
		{
			str = node;
			u_re = dirichlet[2*str];
			u_im = dirichlet[2*str+1];

			di_b[str * 2] = 1.0;
			di_b[str * 2 + 1] = 0.0;

			for (ipls = 0; ipls < npls; ipls++)
			{
				pr[str * 2 + n_nodes_c * ipls * 2] = u_re;
				pr[str * 2 + 1 + n_nodes_c * ipls * 2] = u_im;
			}

			for (k = ig[str]; k <= ig[str + 1] - 1; k++)
			{
				if (is_dof_bound[jg[k]] == false)
				{
					for (ipls = 0; ipls < npls; ipls++)
					{
						pr[jg[k] * 2 + n_nodes_c * ipls * 2] -= ggu_b[k * 2] * u_re - ggu_b[k * 2 + 1] * u_im;
						pr[jg[k] * 2 + 1 + n_nodes_c * ipls * 2] -= ggu_b[k * 2] * u_im + ggu_b[k * 2 + 1] * u_re;
					}
					ggu_b[k * 2] = 0.0;
					ggu_b[k * 2 + 1] = 0.0;
				}
				ggl_b[k * 2] = 0.0;
				ggl_b[k * 2 + 1] = 0.0;
			}
		}
	}

	for (i = 0; i < n_nodes_c; i++)
	{
		for (j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			k = jg[j];
			if (is_dof_bound[k] == true)
			{
				if (is_dof_bound[i] == false)
				{
					for (ipls = 0; ipls < npls; ipls++)
					{
						pr[i * 2 + n_nodes_c * ipls * 2] -= ggl_b[j * 2] * dirichlet[2*k] - ggl_b[j * 2 + 1] * dirichlet[2*k+1];
						pr[i * 2 + 1 + n_nodes_c * ipls * 2] -= ggl_b[j * 2] * dirichlet[2*k+1] + ggl_b[j * 2 + 1] * dirichlet[2*k];
					}
					ggl_b[j * 2] = 0.0;
					ggl_b[j * 2 + 1] = 0.0;
				}
				ggu_b[j * 2] = 0.0;
				ggu_b[j * 2 + 1] = 0.0;
			}
		}
	}

	if(is_dof_bound) {delete [] is_dof_bound; is_dof_bound=NULL;}
	if(dirichlet) {delete [] dirichlet; dirichlet=NULL;}
}
