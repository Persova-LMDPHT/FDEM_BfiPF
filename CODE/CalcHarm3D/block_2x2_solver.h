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
 *  This file contains headers of subroutines for matrix-vector operations with a matrix stored in 2x2-block sparse format          
 *                                                                                                             
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/block_2x2_solver.h                                                 
 *  Novosibirsk State Technical University,                                                                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                      
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                       
 *                                                                                                                                
*/                                                                                                             

#pragma once
#include "base_solver.h"
//------------------------------------------------------------------------

class Block_2x2_solver 
{
public:
	Block_2x2_solver();
	~Block_2x2_solver();

	// multiplication of a single block
	inline void Mult_block_2x2(double *a, int size, double *x, double *y);

	// multiplication of a transposed block 
	inline void Mult_MV_block_2x2_transp(double *a, int size, double *x, double *y);

	// multiplication of a block matrix by a vector
	void Mult_MV_block_2x2(int nb, int *ig, int *jg, int *idi, int *ijg, 
		double *di_block, double *ggl_block, double *x, double *y, double *y_omp);

	// multiplication of a transposed block matrix by a vector
	void Mult_MV_block_2x2_transp(int nb, int *ig, int *jg, int *idi, int *ijg, 
		double *di_block, double *ggl_block, double *x, double *y, double *y_omp);


	///////////////// subroutines for block-diagonal preconditioning ////////////////////

	// Build block-diagonal preconditioner
	int Build_block_diag_preconditioner(int n, int *idi, double *di_block,
		double *df, int *idi_f,
		double *ggl_f, double *ggu_f);

	// Build block-diagonal preconditioner
	int Build_block_diag_preconditioner(int nb, int *idi, 
		double *di_block, double *df, double *ggl_f, double *ggu_f);

	// Build complex-valued block-diagonal preconditioner
	int Build_complex_diag_preconditioner(int nb, int *idi, double *di_block, double *df);

	// Solve lower triangular system with block-diagonal matrix
	int solve_l_blockdiag(int n, double *df, int *idi_f, double *ggl_f, 
		double *f, double *x);

	// Solve upper triangular system with block-diagonal matrix
	int solve_u_blockdiag(int n, double *df, int *idi_f, double *ggl_f, 
		double *f, double *x);

	// Solve lower triangular system with block-diagonal matrix
	int solve_l_blockdiag(int n, double *df, double *ggl_f,double *f, double *x);

	// Solve upper triangular system with block-diagonal matrix
	int solve_u_blockdiag(int n, double *df, double *ggl_f, double *f, double *x);

	// Multiplication of lower triangular system with block-diagonal matrix by a vector
	void mult_l_blockdiag(int nb, double *df, double *ggl_f, double *x, double *y);

	// Multiplication of upper triangular system with block-diagonal matrix by a vector
	void mult_u_blockdiag(int nb, double *df, double *ggu_f, double *x, double *y);

	// Complex-valued incomplete Cholesky factorization
	int LLT_Cmplx(int nb, int *ig, int *jg, int *idi, int *ijg, double *di, double *gg,
		double *d, double *sg);

	// A forward solve with complex-valued sparse triangular matrix
	int SolveL_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
		double *f, double *x);

	// A backward solve with complex-valued sparse triangular matrix
	int SolveU_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
		double *f, double *s, double *x);

	// Permutation of a vector
	void Perm(int nb, double *x, double *y);
};

//------------------------------------------------------------------------