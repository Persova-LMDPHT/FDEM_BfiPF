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
 *  This file contains headers of subroutines for managing the number of OMP threads  
 *                                                                                                                                                   
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/ControlOMP.h                                                                                            
 *  Novosibirsk State Technical University,                                                                                                                                                     
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                                                
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                     
 *  Version 1.3 April 12, 2021                                                                                                                                                          
*/

//------------------------------------------------------------------------
#pragma once
//------------------------------------------------------------------------
class ControlOMP
{
	// One global object of this class is started for the entire program and
	// the management of the number of threads in OpenMP is carried out through it
private:
	bool isInit;
	int NUMBER_OF_THREADS; // current number of threads in OpenMP
	int MAX_THREADS; // maximum possible number of threads in OpenMP
	int nMinDotProduct;
	int nMinSparseMultMV;
	int nMinSparseMultMV2x2;
public:

	void InitNumThreads(); // called once at the beginning of the program
	void SetNumThreads(int num);
	int GetNumberOfThreads();
	int GetMaxThreads();
	int GetNMinDotProduct();
	int GetNMinSparseMultMV();
	int GetNMinSparseMultMV2x2();

	ControlOMP();
	~ControlOMP();
};
