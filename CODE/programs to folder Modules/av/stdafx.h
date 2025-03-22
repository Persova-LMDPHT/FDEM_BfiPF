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
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 * Version 2.0 December 10, 2024
*/

#pragma once

#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <complex>

using namespace std;

#include "For_Solvers.h"

#define PI    3.1415926535897932
#define MU_0  1.2566370614359173e-6

__forceinline istream& operator > (istream& file,double&  data)
{
	file.read((char*)&data,sizeof(data));
	return file;
}

__forceinline ostream& operator < (ostream& file,const double& data)
{
	file.write((char*)&data,sizeof(data));
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

#include "mkl.h"
