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
 *  This file contains code of some basic routines for vector-matrix operations and error messages
 *
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/For_Solvers.cpp         
 *  Modified by D.Sc. Denis V. Vagin                                                                           
 *  Novosibirsk State Technical University,                                                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                          
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                           
 *  Version 2.0 December 10, 2024                                                                              
*/

#include "stdafx.h"
#include "in_out.h"
extern ofstream logfile;
//------------------------------------------------------------------------     
// reports a memory allocation error                                           
//------------------------------------------------------------------------     
void Memory_allocation_error(const char *var, const char *func)
{
	string str;
	str = "MEMORY ALLOCATION ERROR for variable ";
	str = str + '\"' + var + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}

//------------------------------------------------------------------------   
// reports a file open error and throws an exception                         
//------------------------------------------------------------------------   
void Cannot_open_file(const char *fname, const char *func)
{
	string str;
	str = "CANNOT OPEN FILE ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
	cout << str << flush;
	throw logic_error(str);
}
//------------------------------------------------------------------------         
// reports an error opening the file, but the program continues                    
//------------------------------------------------------------------------         
void Cannot_open_file_but_continue(const char *fname, const char *func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
}
//------------------------------------------------------------------------   
// Dot product                                                               
//------------------------------------------------------------------------   

double Scal(double *a, double *b, int n)
{
	int i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}
//------------------------------------------------------------    
// Matrix-vector multiplication in dense format                   
//------------------------------------------------------------    
void Mult_Plot(double *a, double *x, double *y, int n)
{
	int i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*n;
		for(j=0;j<n;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}
//------------------------------------------------------------   
// Euclidian norm of a vector                                    
//------------------------------------------------------------   
double Norm_Euclid(double *a, int n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}
//-------------------------------------------------------------        
// Matrix-vector multiplication in sparse format                       
//-------------------------------------------------------------        
void Mult_MV(int *ig, int *jg, double *ggl, double *ggu, double *di, double *x, double *y, int n)
{
	int i, j, k;

	for(i=0; i<n; i++)
	{
		y[i] = di[i]*x[i];
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			y[i] += ggl[j]*x[k];
			y[k] += ggu[j]*x[i];
		}
	}
}
