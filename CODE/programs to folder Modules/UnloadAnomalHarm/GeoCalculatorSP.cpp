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
 *  This file contains main function for getting points for output primary field
 *  
 *  Written by D.Sc. Denis V. Vagin                                           
 *  Novosibirsk State Technical University,                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                               
 * Version 2.0 December 10, 2024                                                          
*/

#include "stdafx.h"
#include "TaskSP.h"
#include "time_approx.h"

bool CheckStop(void)
{
	bool fstop;
	ifstream ifstop;
	fstop=false;
	ifstop.open("stop");
	if(ifstop){
		fstop=true;
		ifstop.close();
	}
	ifstop.clear();
	return fstop;
}

using namespace std;

ofstream logfile;

extern double Scal(double* a, double* b, int n);

extern bool CheckStop(void);

void Memory_allocation_error(const char* var, const char* func)
{
	string str;
	str = "MEMORY ALLOCATION ERROR for variable ";
	str = str + '\"' + var + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr << str << flush;
	cout << str << flush;
	throw logic_error(str);
}

void Cannot_open_file(const char* fname, const char* func)
{
	string str;
	str = "CANNOT OPEN FILE ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr << str << flush;
	cout << str << flush;
	throw logic_error(str);
}

void Cannot_open_file_but_continue(const char* fname, const char* func)
{
	string str;
	str = "Cannot open file ";
	str = str + '\"' + fname + '\"' + " in function " + '\"' + func + "\"\n";
	logfile << str << flush;
	cerr << str << flush;
}

PointXYZ operator*(const PointXYZ& p, const double& a)
{
	return PointXYZ(p.x * a, p.y * a, p.z * a);
}

PointXYZ operator+(const PointXYZ& p1, const PointXYZ& p2)
{
	return PointXYZ(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
}

ifstream& operator>>(ifstream& inf, PointXYZ& p)
{
	inf >> p.x >> p.y >> p.z;
	return inf;
}

ifstream& operator>(ifstream& inf, PointXYZ& p)
{
	inf > p.x > p.y > p.z;
	return inf;
}

ofstream& operator<<(ofstream& outf, const PointXYZ& p)
{
	outf << p.x << " " << p.y << " " << p.z;
	return outf;
}

ofstream& operator<(ofstream& outf, const PointXYZ& p)
{
	outf < p.x < p.y < p.z;
	return outf;
}

void CopyV(const int& n, const double* from, double* to)
{
	int i;
	for (i = 0; i < n; i++)
		to[i] = from[i];
}

double Scal(const int& n, const double* v1, const double* v2)
{
	double s = 0;
	int i;
	for (i = 0; i < n; i++)
		s += v1[i] * v2[i];
	return s;
}

void Mult_Plot(double* a, double* x, double* y, int n)
{
	int i, j, temp;
	double sum;

	for (i = 0; i < n; i++)
	{
		sum = 0.0;
		temp = i * n;
		for (j = 0; j < n; j++)
			sum += a[temp + j] * x[j];

		y[i] = sum;
	}
}

void AddV(const int& n, const double* from, double* to, const double& mlt = 1)
{
	for (int i = 0; i < n; i++)
		to[i] += mlt * from[i];
}

int CalcSP()
{
	int i;
	ifstream inf;
	PointXYZ pntE0;
	int npntE0, mtrpntE0;
	Time_approx_for_vfem* ta;
	ta = new Time_approx_for_vfem(false, NULL, true, NULL);
	inf.clear();
	inf.open("z_sig_2d");
	if (!inf)
	{
		cout << "Error in open file " << "z_sig_2d" << endl;
		logfile << "Error in open file " << "z_sig_2d" << endl;
		return 1;
	}
	inf >> ta->nzl;
	ta->zlay.resize(ta->nzl);
	for (i = 0; i < ta->nzl; i++)
	{
		inf >> ta->zlay[i];
		inf.ignore(1000, '\n');
	}
	inf.close();
	inf.clear();
	return ta->Read_data();
}

int main(int argc, char **argv)
{
	try
	{
		logfile.open("logUnloadAnomalNodesNonStat");
		if (CalcSP() != 0)
		{
			logfile << "Error!!!" << endl;
			logfile.close();
			return 1;
		}
	}
	catch (...)
	{
		return 1;
	}

	return 0;
}
