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
 *  This file contains main function for calculating the eleptic source part of field for harmonic task
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 * Version 2.0 December 10, 2024
 *                                                                             
*/

#include "stdafx.h"
#include "Harm2D.h"

ofstream logfile;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	try
	{
		logfile.open("LogHarm2d_Ax");

		int i, j, f;
		char tdir[256];
		ofstream outf;
		vector<double> Z, Sigma;
		bool with_gaps=false;
		double freq;

		_mkdir("Ax");

		system("copy mtr3d2d Ax\\mtr3d2d");
		system("copy nthreads.txt Ax\\nthreads.txt");

		ifstream fin;

		fin.open("nu");
		fin >> freq;
		fin.close();
		fin.clear();

		int nfreq=0;
		double curr=1.0;

		fin.open("currentval");
		fin >> curr;
		fin.close();
		fin.clear();

		system("copy inf2tr.dat Ax\\inf2tr.dat");
		system("copy nvtr.dat Ax\\nvtr.dat");
		system("copy nvkat2d.dat Ax\\nvkat2d.dat");
		system("copy rz.dat Ax\\rz.dat");
		system("copy l1.dat Ax\\l1.dat");
		system("copy sigma Ax\\sigma");
		system("copy sigmaZ Ax\\sigmaZ");
		system("copy mu Ax\\mu");
		system("copy dpr Ax\\dpr");
		system("copy mtr3d2d Ax\\mtr3d2d");
		system("copy tsize.dat Ax\\tsize.dat");
		system("copy rz.txt Ax\\rz.txt");
		system("copy r.dat Ax\\r.dat");
		system("copy z.dat Ax\\z.dat");
		system("copy group Ax\\group");
		system("copy sours Ax\\sours");

		_chdir("Ax");

		outf.open("nu");
		outf<<freq<<endl;
		outf.close();
		outf.clear();

		ProcTask2DLine_Harmonic(curr,nfreq);

		_chdir("..");
	}
	catch(...) 
	{ 
		cout << "Unknown exception.\n"  << flush;
		logfile << "Unknown exception.\n"  << flush;
	}

	logfile.close();
	logfile.clear();

	return 0;
}
