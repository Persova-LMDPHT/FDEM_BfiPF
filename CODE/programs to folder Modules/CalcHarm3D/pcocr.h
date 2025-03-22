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
 *  This file contains the header of the COCR solver with folded preconditioner and minimal residual smoothing  
 *  The matix stored in the sparse format                                                                       
 *                                                                                                                 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/pcocr.h                                            
 *  Novosibirsk State Technical University,                                                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                      
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                       
                                          
*/                                                                                                                 
 
#pragma once
//------------------------------------------------------------------------
class PCOCR
{
public:
	int PCOCR_2x2_Folded(int n, int *ig, int *jg, int *idi, int *ijg, double *di, double *gg, double *pr,
		double *x, double eps, int maxiter, double *y_omp,
		int kpar, int n_edges_c, int n_nodes_c,
		double (*xyz)[3], int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
		int *ig_t, int *jg_t, double *gg_t, int *is_node_bound);
	PCOCR();
	~PCOCR();
};
