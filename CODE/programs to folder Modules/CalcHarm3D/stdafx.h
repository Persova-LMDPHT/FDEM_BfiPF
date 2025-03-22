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
 *  This file contains list of include headers and global constants and functions
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024
*/

#pragma once

#include <stdio.h>
#include <tchar.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <stack>
#include <complex>
#include <iomanip>

using namespace std;

#include "in_out.h"
#include "For_Solvers.h"


#define MU_0  1.25663706143591729e-6
#define PI    3.1415926535897932384626433832795
#define DPR_0 8.84194128288307421e-12

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE


__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const float& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,float&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const int& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,int&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const bool& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,bool&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned char& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned char&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const short& data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,short&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const unsigned int & data)
{
	file.write((char*)&data,sizeof(data));
	return file;
}

__forceinline istream& operator > (istream& file,unsigned int &  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

#include "mkl.h"

#include <windows.h>

void reverse_vec3(double *a);
double norma3(double *a);
void normalize3(double *a);
double mult_scal3(double *a,double *b);
void mult_vec3(double *a,double *b,double *c);
void mult_matrix33(double a[][3],double b[][3],double c[][3]);
double GetDeterminant33(double m[][3]);
void TransposeMatrix33(double a[][3]);
void GetRotationMatrix33(double M_T[][3],double *v,double angle);
void InverseMatrix33(double m[][3],double m_1[][3]);
void InverseMatrix33WithDet(double m[][3],double m_1[][3],double det);
