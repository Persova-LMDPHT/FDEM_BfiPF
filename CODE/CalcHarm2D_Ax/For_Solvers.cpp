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
#include "For_Solvers.h"
#include "in_out.h"
extern ofstream logfile;

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

void Cannot_open_file_but_continue(const char *fname, const char *func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr    << str << flush;
}

double Scal(double *a, double *b, long n)
{
	long i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}

void Mult_Plot_AV(double *a, double *x, double *y, long n, long m)
{
	long i, j, temp;
	double sum;

	for(i=0;i<n;i++)
	{
		sum = 0.0;
		temp = i*m;
		for(j=0;j<m;j++)
			sum += a[temp+j]*x[j];

		y[i] = sum;
	}
}

double Projection_On_Axis(double *v,double *o)
{
	double value;

	value = Scal(v,o,3)/sqrt(Scal(o,o,3));

	return value;
}

void Mult_Plot(double *a, double *x, double *y, long n)
{
	long i, j, temp;
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

double Norm_Euclid(double *a, long n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}

double Norm_Max(double *a, long n)
{
	double max, current;
	long i;

	max = 0.0;

	for(i=0; i<n; i++)
	{
		current = a[i]*a[i];
		if(current > max)
			max = current;
	}

	return max;
}

void Mult_MV(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n)
{
	long i, j, k;

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

double Spline(double x, long n, double *xyz, double *values)
{
	double s, xi;
	long i, t, flag;

	flag = 0;

	if (x < xyz[0])
	{
		return values[0];
	}

	if (x > xyz[n-1])
	{
		return values[n-1];
	}

	for(i=0; i<n-1; i++)
	{
		if(x >= xyz[i]  &&  x <= xyz[i+1])
		{
			t = i;
			flag = 1;
			break;
		}
	}

	if(flag == 1)
	{
		xi = (x - xyz[t])/(xyz[t+1] - xyz[t]);
		s = (1.0 - xi)*values[t] + xi*values[t+1];
	}
	else
	{
		s = 0.0;
	}

	return s;
}
