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
 *  This file contains headers of subroutines for basic vector and matrix-vector operations
 *                                                                                         
 *  Written by D.Sc. Denis V. Vagin, Prof. Marina G. Persova                    
 *  Novosibirsk State Technical University,                                     
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                           
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)            
 * Version 2.0 December 10, 2024                                                
 *                                                                      
 */                                


#pragma once

class Base_solver
{
// The Base_solver class contains elementary routines for constructing iterative solvers (dot product, vector norm, etc.)
public:
	Base_solver();
	~Base_solver();

	inline double Scal(double *x, double *y, int n);   // dot product
	double Norm_Euclid(double *x, int n);              // vector euclidean norm
	double Projection(double *vec, double *axis);      // projection of a vector onto an axis
	double Relative_Error(double *analytic, double *numeric, int n); // relative error
	double Spline(double x, int n, double *xyz, double *values);	 // linear interpolation

	// Matrix-vector multiplication  (in dense format)
	void Mult_Plot(double *a,double *pr,double *rez,int n); 

	// Givens rotation
	int  Givens1(double& x, double& y, double& c, double& s);
	void Givens2(double& x, double& y, double c, double s);
	int  Givens(double *a, double *f, int n);

	// Solution of a system with a lower triangular matrix in a dense format
	int Undirect(double *a, double *b, double *x, int n);

	// Solution of a SLAE with a square matrix whose lower triangle contains only one non - zero subdiagonal.	
	// This is necessary if the Arnoldi orthogonalization fails
	int Solve_square_subdiag(double *a, double *b, double *x, int n);

	// writing to the file the relative discrepancy with which they came out, eps, the number of iterations and the time of solving the SLAE
	int WriteKitChrono(char *fname, double residual, double eps, int iter, double time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, int n, int size_jg);
	int Write_kit(char *fname, double residual, double eps, int iter, int time, double change_of_solution);

	// for complex-valued arithmetic (vectors are stored in the usual format, complex numbers (scalars) - in std::complex<double>)

	// dot product for complex-valued vectors
	std::complex<double> ScalCmplxTrue(double *x, double *y, int nb); 

	// complex conjugate dot product
	std::complex<double> ScalCmplx(double *x, double *y, int nb); 

	// multiplication of a vector by a complex number
	void MultCmplxNumVect(std::complex<double> a, double *x, double *y, int nb);

	// multiplication of components of one complex vector by components of another complex vector
	void MultCmplxVectVect(int nb, double *a, double *b, double *c);

	// division of the components of one complex vector into components of another complex vector
	void DivCmplxVectVect(int nb, double *a, double *b, double *c);

	// z = x + a*y
	void Cmplx_axpy(std::complex<double> a, double *x, double *y, double *z, int nb);

	double GetMinimiz(int n, double *ap, double *r);
};

//----------------------------------------------------
