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
 *  This file contains code for using MKL PARDISO solver
 *  
 *  
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/pardiso.cpp  
 *  Modified by D.Sc. Denis V. Vagin                                                                
 *  Novosibirsk State Technical University,                                                         
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                               
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                
 *  Version 2.0 December 10, 2024                                                                   
*/

#include "stdafx.h"
#include "pardiso.h"

extern ofstream logfile;
//-----------------------------------------------------------     
// Constructor                                                    
//-----------------------------------------------------------     
pardiso_solver::pardiso_solver()
{
	ia=NULL;
	ja=NULL;
	a=NULL;
	perm=NULL;

	n = 0;

	init();
}
//-----------------------------------------------------------   
// Destructor                                                   
//-----------------------------------------------------------   
pardiso_solver::~pardiso_solver()
{
	clear();
}
//-----------------------------------------------------------  
// Pardiso parameters initialization for complex matrix        
//-----------------------------------------------------------  
void pardiso_solver::init()
{
	int i;

	mtype = 13; //2;
	nrhs = 1;
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	phase = 13;

	for(i=0;i<64;i++){pt[i]=0;}
	for(i=0;i<64;i++){iparm[i]=0;}

	msglvl=0;
}

void pardiso_solver::clear()
{
	if(a){delete [] a;  a=NULL;}
	if(ia){delete [] ia; ia=NULL;}
	if(ja){delete [] ja; ja=NULL;}
	if(perm){delete [] perm; perm=NULL;}
}
//-----------------------------------------------------------       
// Converting and factorization for complex matrix                  
//-----------------------------------------------------------       
void pardiso_solver::factorize(int nb,int *ig,int *jg,double *ggl,double *ggu,double *di,int *idi,int *ijg,int nthreads)
{
	int i;
	int ig_n_1=ig[nb];
	int sz_iptr=0;
	int sz_jptr=0;

	init();

	if (n == 0)
	{
		clear();
	}

	mkl_set_num_threads(nthreads);

	f.From2x2ToCSR_Complex_1_NS(nb, ig, idi, ijg, &sz_iptr, &sz_jptr);

	if (ia==NULL) ia = new MKL_INT64[sz_iptr];
	if (ja==NULL) ja = new MKL_INT64[sz_jptr];
	if (a==NULL)   a = new double[sz_jptr*2];

	f.From2x2ToCSRComplex_2_NS(nb, ig, jg, idi, ijg, di, ggl, ggu, ia, ja, a);

	for(i=0;i<sz_iptr;i++){ia[i]++;}
	for(i=0;i<sz_jptr;i++){ja[i]++;}

	if (perm==NULL)	perm = new MKL_INT64[nb];


	phase = 12;
	n = nb;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
	cout << "factorize info=" << info << '\n';
}
//-----------------------------------------------------------  
// Forward reverse                                             
//-----------------------------------------------------------  
void pardiso_solver::solve_nrhs(int _nrhs,double *pr,double *x)
{
	int i;
	for(i=0; i<n*2; i++) {x[i]=0.0;}
	phase = 33;
	nrhs = _nrhs;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, pr, x, &info);
}
//----------------------------------------------------------- 
// Finish PARDISO                                             
//----------------------------------------------------------- 
void pardiso_solver::stop_solver()
{
	phase = -1;
	pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, NULL, NULL, &info);
}
