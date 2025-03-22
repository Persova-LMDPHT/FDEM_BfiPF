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
 *  This file contains code for calculating 2D tasks.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024    
*/

#include <io.h>
#include <fcntl.h>
#include <share.h>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

char str[1024];

ofstream logfile;

void close_logfile()
{
	logfile.close();
	logfile.clear();
}

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (!(retp=CreateProcessA(NULL,(LPSTR)(const char*)cmdline,NULL,NULL,FALSE,0,NULL,workdir,&si,&pi)))
	{
		sprintf(str,"Can't create process for %s, error code %d",cmdline,GetLastError());
		cout<<str<<endl;
		logfile<<str<<endl;
		return 1;
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
	GetExitCodeProcess(pi.hProcess, (LPDWORD)&retp);
	CloseHandle(pi.hProcess);
	CloseHandle(pi.hThread);
	return retp;
}

void run(char *cmd)
{
	int retp;
	retp=CreateProcessForEXE(cmd,NULL);
	if(retp)
	{
		cout<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		close_logfile();
		exit(retp);
	}
}

struct PointXYZ
{
	double x, y, z;
	PointXYZ()
	{
		x = y = z = 0;
	}
	PointXYZ(const double& _x, const double& _y, const double& _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	const PointXYZ& operator+=(const PointXYZ& p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}
	const PointXYZ& operator-=(const PointXYZ& p)
	{
		x -= p.x;
		y -= p.y;
		z -= p.z;
		return *this;
	}
	void read(ifstream& inf)
	{
		inf >> x >> y >> z;
	}
	void write(ofstream& ofp)
	{
		ofp << x << ' ' << y << ' ' << z;
	}
};

struct LineXYZ
{
	PointXYZ A, B;
	void read(ifstream& inf)
	{
		A.read(inf);
		B.read(inf);
	}
	void write(ofstream& ofp)
	{
		A.write(ofp);
		ofp << ' ';
		B.write(ofp);
		ofp << '\n';
	}
};

int main(int argc,char **argv)
{
	int i,j,k,retp,ipls,npls,nprof,p1,p2,met2d;
	ifstream inf;
	ofstream ofp;
	int fmet2d;
	vector<LineXYZ> vGels;

	logfile.open("logfile2D");

	npls = 0;
	inf.open("group");
	if (!inf)
	{
		logfile << "Error in open file " << "group" << endl;
		cout << "Error in open file " << "group" << endl;
		return 1;
	}
	inf >> nprof;
	for (i = 0; i < nprof; i++)
	{
		inf >> k >> p1 >> p2;
		npls += p2 - p1 + 1;
	}
	inf.close();
	inf.clear();

	inf.open("sours");
	if (!inf)
	{
		logfile << "Error in open file " << "sours" << endl;
		cout << "Error in open file " << "sours" << endl;
		return 1;
	}
	vGels.resize(npls);
	for (i = 0; i < npls; i++)
	{
		vGels[i].read(inf);
	}
	inf.close();
	inf.clear();

	ofp.open("srsclca");
	for (i = 0; i < npls; i++)
	{
		vGels[i].write(ofp);
	}
	ofp.close();
	ofp.clear();

	ofp.open("clcnplsa");
	ofp << npls << '\n';
	ofp.close();
	ofp.clear();

	ofp.open("srsclcgsz");
	for (i = 0; i < npls; i++)
	{
		ofp << 1 << '\n';
	}
	ofp.close();
	ofp.clear();

	ofp.open("srsvala");
	for (i = 0; i < npls; i++)
	{
		ofp << 1.0 << '\n';
	}
	ofp.close();
	ofp.clear();

	logfile<<"start Meshing"<<endl;
	run("..\\..\\Modules\\Mesh2D_FD.exe");
	logfile<<"finish Meshing"<<endl;

	met2d = 0;
	inf.open("met2d");
	if (inf)
	{
		inf >> met2d;
		inf.close();
	}
	inf.clear();

	if (!met2d) 
	{
		logfile << "start Solving 2d" << endl;
		run("..\\..\\Modules\\CalcHarm2DHEL_U.exe");
		logfile << "finish Solving 2d" << endl;
	}
	else
	{
		system("del /q inf2tr.dat nvtr.dat nvkat2d.dat rz.dat l1.dat tsize.dat rz.txt r.dat z.dat currentval");

		system("copy inf2tr.dat1 inf2tr.dat");
		system("copy nvtr.dat1 nvtr.dat");
		system("copy nvkat2d.dat1 nvkat2d.dat");
		system("copy rz.dat1 rz.dat");
		system("copy l1.dat1 l1.dat");
		system("copy tsize.dat1 tsize.dat");
		system("copy rz.txt1 rz.txt");
		system("copy r.dat1 r.dat");
		system("copy z.dat1 z.dat");
		system("copy currentval1 currentval");

		logfile << "start Solving 2d" << endl;
		run("..\\..\\Modules\\CalcHarm2DHEL_AV.exe");
		logfile << "finish Solving 2d" << endl;
	}

	close_logfile();

	return 0;
}
