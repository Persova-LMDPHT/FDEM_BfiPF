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
 *  This file contains headers of subroutines for the folded preconditioner                     
 * 
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/FoldedPreconditioner.h                                                                                       
 *  Modified by D.Sc. Denis V. Vagin                                                  
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024                                                                                                          
*/                                                                                        

#pragma once
//------------------------------------------------------------------------
class FoldedPreconditioner
{
public:
	int n; // size of the original system
	int m; // space dimension

	double zero;

	int *iptr_g;
	int *jptr_g;
	double *gg_g;

	int *iptr_gt;
	int *jptr_gt;
	double *gg_gt;

	double *diag_a;
	double *diag_gagt;
	double *help;
	double *help2;
	double *help3;

	FoldedPreconditioner();
	~FoldedPreconditioner();

	// Construction of the discrete gradient matrix
	void BuildGMatrix(int kpar, int n_edges_c, int n_nodes_c,
		int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound);

	// Initialization of the folded preconditioner
	void Prepare(int kpar, int n_edges_c, int n_nodes_c,
		int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound,
		int *ig, int *jg, double *di, double *gg);

	// Build the main diagonal of the matrix P(G+M)PT
	void BuildDiag(int *ig, int *jg, double *di, double *gg);

	// Applying the folded preconditioner to a vector
	void ApplyPreconditioner(double *x, double *y);

	// Build the main complex-valued diagonal of the matrix P(G+iwM)PT
	void BuildDiagComplex(int *ig, int *jg, int *idi, int *ijg, double *di, double *gg);

	// Applying the complex-valued folded preconditioner to a vector
	void ApplyPreconditionerComplex(double *x, double *y);

	// Conversion to the CSR format
	void FromRSFtoCSR_1(int n, int *ig, int *sz_iptr, int *sz_jptr);

	// Conversion to the CSR format
	void FromRSFtoCSR_2(int n, int *ig, int *jg, double *di, double *ggl, double *ggu,
		int *iptr, int *jptr, double *aelem);

	// Conversion to the block CSR format
	void From2x2ToCSR2x2_1(int n, int *ig, int *idi, int *ijg,
		int *sz_iptr, int *sz_jptr, int *sz_ijptr, int *sz_aelem);

	// Conversion to the block CSR format
	void From2x2ToCSR2x2_2(int n, int *ig, int *jg, int *idi, int *ijg,
		double *di_block, double *gg_block, 
		int *iptr, int *jptr, int *ijptr, double *aelem);
};
