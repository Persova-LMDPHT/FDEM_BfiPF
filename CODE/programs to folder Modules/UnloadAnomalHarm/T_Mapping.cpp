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
 *  This file contains code for working with edge mesh in 3D VFEM
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                           
 *  Novosibirsk State Technical University,                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                         
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)          
 * Version 2.0 December 10, 2024                                              
*/

#include "stdafx.h"
#include "T_Mapping.h"
#include "in_out.h"

extern ofstream logfile;
extern void Memory_allocation_error(const char *var, const char *func);
extern void Cannot_open_file(const char *fname, const char *func);

double Scal(double *a, double *b, int n)
{
	long i;
	double sum = 0.0;

	for(i=0; i<n; i++)
		sum += a[i]*b[i];

	return sum;
}

double Norm_Euclid(double *a, int n)
{
	double value;

	value = sqrt(Scal(a,a,n));

	return value;	
}

double Interval(double *x, double *y)
{
	double r;

	r = sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2]));

	return r;
}

double Projection_On_Axis(double *v,double *o)
{
	double value;

	value = Scal(v,o,3)/sqrt(Scal(o,o,3));

	return value;
}

T_Mapping_Vec::T_Mapping_Vec()
{
	this->nver = 0;
	this->xyz = 0;
	this->kuzlov = 0;
	this->kpar = 0;

	ed = NULL;
	edges = NULL;
}
T_Mapping_Vec::T_Mapping_Vec(int (*nver)[14], double (*xyz)[3], int kuzlov, int kpar)
{
	this->nver = nver;
	this->xyz = xyz;
	this->kuzlov = kuzlov;
	this->kpar = kpar;

	ed = NULL;
	edges = NULL;
}

T_Mapping_Vec::~T_Mapping_Vec()
{
	if(ed) {delete [] ed; ed = NULL;}
	if(edges) {delete [] edges; edges = NULL;}
}
