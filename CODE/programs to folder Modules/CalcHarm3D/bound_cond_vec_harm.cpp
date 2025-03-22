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
 *  The file contains code of functions of the class "Bound_cond_vec_harm"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/bound_cond_vec_harm.cpp
 *  Modified by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024                                                              
*/

#include "stdafx.h"
#include "bound_cond_vec_harm.h"
//-----------------------------------------------------------------------------         
// Constructor for the class Bound_cond_vec_harm (single potentials A-V)                
//-----------------------------------------------------------------------------         
Bound_cond_vec_harm::Bound_cond_vec_harm(long n_nodes, long n_edges, long n_bound_nodes,
						   long *bound_nodes,  long (*edges)[2], long n_elem)
{
	this->n_nodes = n_nodes;
	this->n_edges = n_edges;
	this->n_bound_nodes = n_bound_nodes;
	this->bound_nodes = bound_nodes;
	this->edges = edges;
	this->n_elem = n_elem;

	isBlockBound = m_bound = NULL;

	MakeListBoundsBlocks();
}
//-----------------------------------------------------------------------------  
// Desstructor for the class Bound_cond_vec_harm                                 
//-----------------------------------------------------------------------------  
Bound_cond_vec_harm::~Bound_cond_vec_harm()
{
	if(isBlockBound) {delete [] isBlockBound; isBlockBound=NULL;}
	if(m_bound) {delete [] m_bound; m_bound=NULL;}
}
//------------------------------------------------------------------------    
// Setting the flags if the boundary condition is imposed to nodes and edges  
//------------------------------------------------------------------------    
void Bound_cond_vec_harm::MakeListBoundsBlocks()
{
	long i;
	bool *is_node_bound=NULL;

	if((is_node_bound = new bool[n_nodes])==0) 
		Memory_allocation_error("is_node_bound","Bound_cond_vec_harm::Make_list_of_bound_edges_harm");

	for(i=0; i<n_nodes; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = true;

	if((isBlockBound = new bool[n_edges])==0)
		Memory_allocation_error("isBlockBound","Bound_cond_vec_harm::Make_list_of_bound_edges_harm");;

	for(i=0; i<n_edges; i++)
		isBlockBound[i] = false;

	for(i=0; i<n_edges; i++)
		if(is_node_bound[edges[i][0]]==1 && is_node_bound[edges[i][1]]==1)		
			isBlockBound[i] = true;

	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
//---------------------------------------------------------------------
// Taking into account the first boundary conditions 
//---------------------------------------------------------------------
void Bound_cond_vec_harm::SetHomogenDirichletCond(long *ig, long *jg, long *idi, long *ijg,  
												  double *di, double *gg, double *pr,int npls)
{
	long i, j, k, adr, ipls;
	long size;

	cout << "Set_dirichlet_cond_harm...\n";

	for (i=0; i<n_edges; i++)
	{
		if (isBlockBound[i])
		{
			di[idi[i]] = 1.0;
			size = idi[i+1] - idi[i];
			if(size==2)
				di[idi[i]+1] = 0.0;

			for (j=ig[i]; j<=ig[i+1]-1; j++)
			{
				k = jg[j];
				adr = ijg[j];
				gg[adr] = 0.0;
				size = ijg[j+1] - adr;
				if(size==2)
					gg[adr+1] = 0.0;
			}

			for(ipls=0;ipls<npls;ipls++)
			{
				pr[i*2+n_edges*2*ipls] = pr[i*2+1+n_edges*2*ipls] = 0.0;
			}
		}

		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if (isBlockBound[k])
			{
				adr = ijg[j];
				gg[adr] = 0.0;
				size = ijg[j+1] - adr;
				if(size==2)
					gg[adr+1] = 0.0;
			}
		}
	}
}
//---------------------------------------------------------------------
// Taking into account the first boundary conditions only for the right-hand side 
//---------------------------------------------------------------------
void Bound_cond_vec_harm::SetHomogenDirichletCondPr(long *ig, long *jg, long *idi, long *ijg,  
												  double *di, double *gg, double *pr)
{
	long i, j, k, adr;
	long size;

	for (i=0; i<n_edges; i++)
	{
		if (isBlockBound[i])
		{
			pr[i*2] = pr[i*2+1] = 0.0;
		}
	}
}
