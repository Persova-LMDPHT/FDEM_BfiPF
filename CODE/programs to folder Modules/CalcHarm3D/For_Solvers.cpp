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
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/For_Solvers.cpp
 *  Modified by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024
 *  
*/

#include "stdafx.h"
#include "T_Brick.h"
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
//------------------------------------------------------------ 
// Projection of the vector v onto the o axis                  
//------------------------------------------------------------ 
double Projection_On_Axis(double *v,double *o)
{
	double value;

	value = Scal(v,o,3)/sqrt(Scal(o,o,3));

	return value;
}
//------------------------------------------------------------ 
// Matrix-vector multiplication in dense format                
//------------------------------------------------------------ 
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
//------------------------------------------------------------  
// Euclidian norm of a vector                                   
//------------------------------------------------------------  
double Norm_Euclid(double *a, long n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}
//------------------------------------------------------------
// Max norm of a vector                                       
//------------------------------------------------------------
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
//-------------------------------------------------------------    
// Matrix-vector multiplication in sparse format                   
//-------------------------------------------------------------    
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
//-----------------------------------------------------------       
// Relative error                                                   
//-----------------------------------------------------------       
double Relative_Error(double *analytic, double *numeric, long n)
{
	double *razn=NULL;
	double norm_analytic;
	double norm_razn;
	double error;
	long i;

	razn = new double[n];
	if(razn == 0)
		Memory_allocation_error("razn", "Relative_Error");

	for(i=0;i<n;i++)
		razn[i] = analytic[i] - numeric[i];

	In_Out R;
	R.Write_Txt_File_Of_Double("razn", razn, n, 1);

	norm_analytic = Norm_Euclid(analytic, n);
	norm_razn = Norm_Euclid(razn, n);

	error = norm_razn/norm_analytic;

	if(razn) {delete [] razn; razn=NULL;}

	return error;
}
//-----------------------------------------------------------  
// Maximum of 2 numbers                                        
//-----------------------------------------------------------  
long Max_Long(long a, long b)
{
	if(a > b) return a;
	return b;
}
//-----------------------------------------------------------     
// Minimum of 2 numbers                                           
//-----------------------------------------------------------     
long Min_Long(long a, long b)
{
	if(a < b) return a;
	return b;
}
//-----------------------------------------------------------   
// Sort 2 numbers                                               
//-----------------------------------------------------------   
void Sort2(long *a, long *b)
{
	long tmp;
	if(*a < *b)
	{
		tmp = *a;
		*a = *b;
		*b = tmp;
	}
}
//----------------------------------------------------------- 
// Distance between two 3D-points                             
//----------------------------------------------------------- 
double Interval(double *x, double *y)
{
	double r;

	r = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]));

	return r;
}
//----------------------------------------------------------- 
// Distance between two parallel lines in 3D                  
//----------------------------------------------------------- 
double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1)
{
	double ax, ay, az;  // direction vector for the 1st line
	double x1, y1, z1;  // point on the first line 
	double x0, y0, z0;  // point on the second line
	double det1, det2, det3;
	double answer;

	ax = a1[0] - a0[0];
	ay = a1[1] - a0[1];
	az = a1[2] - a0[2];

	x1 = a0[0];
	y1 = a0[1];
	z1 = a0[2];

	x0 = b0[0];
	y0 = b0[1];
	z0 = b0[2];

	det1 = ay*(z1 - z0) - az*(y1 - y0);
	det2 = az*(x1 - x0) - ax*(z1 - z0);
	det3 = ax*(y1 - y0) - ay*(x1 - x0);

	answer = sqrt(det1*det1 + det2*det2 + det3*det3)/sqrt(ax*ax + ay*ay + az*az);

	return answer;
}
//-----------------------------------------------------------      
// Linear interpolation                                            
//-----------------------------------------------------------      
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
double Calc_dof(double *J, double *func, long n_local_edge)
{
	double Jt[3];
	Mult_Plot(J, (double*)TANGENT_VECTORS_ON_REFERENCE_CUBE[n_local_edge], Jt, 3);
	return Scal(func, Jt, 3);
}
