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
 *  The reverse communication interface for the iterative solvers                         
 *                                                                                        
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/rci.h                              
 *  Novosibirsk State Technical University,                                                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                  
*/                                                                                        

#pragma once
#include "base_solver.h"
//------------------------------------------------ -----------------------
#define REQ_OK 0 // successful completion of the solver

#define REQ_MULT_MV 1 // need to multiply a matrix by a vector
#define REQ_MULT_L 2 // need to multiply the vector by the left preconditioner
#define REQ_MULT_U 3 // need to multiply the vector by the right preconditioner
#define REQ_SOLVE_L 4 // need to solve system with left preconditioner
#define REQ_SOLVE_U 5 // need to solve system with right preconditioner
#define REQ_PRECOND 6 // need to solve a preconditioned system
#define REQ_STOP_TEST 7 // need to run the solver stop criterion
#define REQ_X0_TEST 8 // need to check if the initial approximation is already a solution

#define REQ_MULT_MV_TRANSP 9 // need to multiply transposed matrix by vector
#define REQ_SOLVE_L_TRANSP 10 //
#define REQ_SOLVE_U_TRANSP 11 //

#define REQ_GET_VECTOR 12 // you need to give the next vector (direction to minimize)

#define REQ_X0_NULL -1 // initial guess is already a solution
#define REQ_DIV_BY_ZERO -2 // division by zero happened inside the solver
#define REQ_MAXITER -3 // emergency exit when the maximum number of iterations is reached

//------------------------------------------------ -----------------------
class RCI
{
protected:
	// variables required for the functioning of the RCI mechanism

	int request; // parameter informing about the result of the solver:
	// request<0 - the solver exited with an error;
	// request=REQ_OK - successful completion of work.
	// A positive value indicates that the user must perform certain actions.

	int stage; // contains information about what stage the solver is at

	double** in; // vector by which to multiply the matrix or act as a preconditioner
	double** out; // result of multiplying a matrix by a vector or a preconditioned vector

public:

	// variables that are in all solvers
	int n; // SLAE dimension
	int nb; // block dimension of SLAE
	double* x; // approximation of the solution vector at the current iteration
	double* pr; // right side vector
	Base_solver bs;


	// variables needed for solver stop criteria
	int iter; // current iteration number
	int maxiter; // maximum number of iterations
	double eps; // SLAE solution accuracy (relative discrepancy norm)
	double r_old; // norm of the true discrepancy before the start of counting
	double eps_x0; // precision to check if the initial guess is already a solution
	double eps_zero; // compare with zero for emergency exit
	double residual; // residual norm at the current iteration
	double residualRel; // norm of relative discrepancy
	double eps_gmres;

public:
	virtual int Run() = 0;

	void PrintIterResidual();
	void DoStopTest(double* r, int* req);
	void DoStopTest(double r, int* req);
	void DoX0Test(double* r, int* req);

	RCI(int n, int maxiter, double eps, double* x, double* pr, double** in, double** out);
	~RCI();
};
