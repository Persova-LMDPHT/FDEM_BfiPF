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
 *  The file contains headers for the functions of the class "Bound_cond_vec_harm"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,                                                                         
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                               
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                
 *  Version 2.0 December 10, 2024                                                                                   
 *  
*/

#pragma once

class Bound_cond_vec_harm
{
public:
	bool *isBlockBound;
	long *bound_nodes;
	long (*edges)[2];

	long n_edges;  
	long n_nodes; 
	long n_elem;
	long n_bound_nodes;

	long n_bound_edges;

	Bound_cond_vec_harm(long n_nodes, long n_edges, long n_bound_nodes,
		long *bound_nodes,  long (*edges)[2], long n_elem);
	~Bound_cond_vec_harm();

	void MakeListBoundsBlocks();
	// Taking into account the first boundary conditions
	void SetHomogenDirichletCond(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr, int npls);
	// Taking into account the first boundary conditions only for the right-hand side
	void SetHomogenDirichletCondPr(long *ig, long *jg, long *idi, long *ijg,  
		double *di, double *gg, double *pr);

	bool *m_bound;
	long *nded;
	long *nded_type;
	long *nvkat;
	long (*nver)[14];
	long (*ed)[25];
	long *nodes_position_in_nded;
	long *edges_position_in_nded;
	long n_nodes_c;
	long n_edges_c;
	long unk_c;
};
