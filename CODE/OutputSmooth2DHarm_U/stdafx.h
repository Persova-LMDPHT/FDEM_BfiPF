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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<windows.h>
#include<iostream>
#include<fstream>
#include<iomanip>
#include<vector>
#include<omp.h>

#define RETCODE_OK				0x0000
#define RETCODE_NOMEM			0x0001
#define RETCODE_NOFILE			0x0002
#define RETCODE_OUTOFRANGE		0x0004
#define RETCODE_SQFROMNEG		0x0008
#define RETCODE_DEVBYZERO		0x0010
#define RETCODE_NOTINIT			0x0020
#define RETCODE_BADFILE			0x0040
#define RETCODE_ERROR			0x0080
#define RETCODE_NOANOMALOBJECTS	0x0100

using namespace std;

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

int sig(const double &v);

int FindInDoubleVec(vector<double> &vec,double elem);
int FindIntervalInDoubleMas(double *vec,int size,double elem);
int FindValueInDoubleMas(double *vec,int size,double elem);

#define PI 3.1415926535897932
#define MU0 1.25663706143592e-6
