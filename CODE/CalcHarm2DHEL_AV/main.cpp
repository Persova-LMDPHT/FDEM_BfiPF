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

#include <sys/stat.h>
#include <io.h>
#include <fcntl.h>
#include <share.h>
#include <windows.h>
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <direct.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

ofstream logfile;

bool isFileExists(char *fname)
{
	ifstream inf;
	bool flag;
	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();
	return flag;
}

vector<double> vin1, vin2, vin3, vout;

struct PointXYZ
{
	double x, y, z;
	PointXYZ() 
	{ 
		x=y=z=0; 
	}
	PointXYZ(const double& _x, const double& _y, const double& _z) 
	{
		x=_x;
		y=_y;
		z=_z; 
	}
	const PointXYZ& operator+=(const PointXYZ& p)
	{
		x+=p.x;
		y+=p.y;
		z+=p.z;
		return *this;
	}
	const PointXYZ& operator-=(const PointXYZ& p)
	{
		x-=p.x;
		y-=p.y;
		z-=p.z;
		return *this;
	}
};

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	char str[256];
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (!(retp=CreateProcessA(NULL,(LPSTR)(const char*)cmdline,NULL,NULL,FALSE,0,NULL,workdir,&si,&pi)))
	{
		sprintf(str,"Can't create process for %s, error code %d",cmdline,GetLastError());
		logfile<<str<<endl;
		return 1;
	}
	WaitForSingleObject(pi.hProcess,INFINITE);
	GetExitCodeProcess(pi.hProcess,(LPDWORD)&retp);
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
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile.close();
		exit(retp);
	}
}

int AddDoubleFile(char *FileName0,char *FileName1,char *FileName2,char *FileName3,int nPntE0,vector<PointXYZ> &TgCompE0,int npls)
{
	int i,j,n,ipls;
	const int size_d=sizeof(double);
	const int size_f=sizeof(float);
	double f1,f2,f3,fs;
	float tmpf;
	FILE *fpi1,*fpi2,*fpi3,*fpo;
	PointXYZ E0;

	n = 3 * nPntE0;

	if (vout.size() < n)
	{
		vin1.resize(n);
		vin2.resize(n);
		vin3.resize(n);
		vout.resize(n);
	}

	fpi1=fopen(FileName0,"rb");
	if(!fpi1)return 1;
	fpi2 = fopen(FileName1, "rb");
	if (!fpi2)return 1;
	fpi3 = fopen(FileName2, "rb");
	if (!fpi3)return 1;
	fpo = fopen(FileName3, "wb");

	for(ipls=0;ipls<npls;ipls++)
	{
		fread(&(vin1.front()), size_d, n, fpi1);
		fread(&(vin2.front()), size_d, n, fpi2);
		fread(&(vin3.front()), size_d, n, fpi3);

		for(i=0;i<n;i++)
		{
			vout[i] = vin1[i] + vin2[i] + vin3[i];
		}

		for (i = 0; i < nPntE0; i++)
		{
			j = 3 * i;
			E0.x = vout[j];
			E0.y = vout[j+1];
			E0.z = vout[j+2];
			tmpf = (float)(E0.x * TgCompE0[i].x + E0.y * TgCompE0[i].y + E0.z * TgCompE0[i].z);
			fwrite(&(tmpf), size_f, 1, fpo);
		}
	}

	fclose(fpi1);
	fclose(fpi2);
	fclose(fpi3);
	fclose(fpo);

	return 0;
}

int SumHarmResFile(char *FileName0,char *FileName1,char *FileName2,char *FileName3)
{
	ifstream inf1,inf2,inf3;
	ofstream ofp;
	int i;
	double t1,t2,t3;
	inf1.open(FileName0);
	if(!inf1)return 1;
	inf2.open(FileName1);
	if(!inf2)return 1;
	inf3.open(FileName2);
	if(!inf3)return 1;
	ofp.open(FileName3);
	ofp<<scientific<<setprecision(14);
	i=0;
	while(!inf1.eof() && !inf2.eof() && !inf3.eof())
	{
		inf1>>t1;
		if(inf1.eof() || !inf1.good())break;
		inf2>>t2;
		if(inf2.eof() || !inf2.good())break;
		inf3>>t3;
		if(inf3.eof() || !inf3.good())break;
		ofp<<(t1+t2+t3);
		i++;
		if(i%6)
		{
			ofp<<' ';
		}
		else
		{
			ofp<<'\n';
		}
	}
	if(i%6)return 1;
	inf1.close();
	inf1.clear();
	inf2.close();
	inf2.clear();
	inf3.close();
	inf3.clear();
	ofp.close();
	ofp.clear();
	return 0;
}

int main(int argc, char **argv)
{
	char buf[1024];
	int i,nPntE0,ipls,npls;
	ifstream inf;
	ofstream ofp;
	vector<int> RecvPlsSzB, RecvPlsSzE;
	double FldHrm[6];
	
	logfile.open("log2dHarmAV");

	logfile<<"start Ax"<<endl;

	run("..\\..\\Modules\\CalcHarm2D_Ax.exe");
	run("..\\..\\Modules\\OutputSmoothAV2DHarm_Ax.exe");

	logfile<<"finish Ax"<<endl;

	if(!isFileExists("Ax\\currentval"))
	{
		system("copy currentval Ax");
	}

	logfile<<"start AV"<<endl;
	

	run("..\\..\\Modules\\av.exe");
	run("..\\..\\Modules\\OutputSmoothAV2DHarm_Er.exe");
	run("..\\..\\Modules\\OutputSmoothAV2DHarm_Ez.exe");

	logfile<<"finish VEL"<<endl;

	npls = 0;

	inf.open("clcnplsa");
	if (!inf)
	{
		logfile << "Error in open file " << "clcnplsa" << endl;
		cout << "Error in open file " << "clcnplsa" << endl;
		return 1;
	}
	inf >> npls;
	inf.close();
	inf.clear();

	if(SumHarmResFile("e2d_ax","e2d_er","e2d_ez","e2d"))
	{
		logfile<<"Error while summing "<<"e2d"<<endl;
		exit(1);
	}

	if(SumHarmResFile("b2d_ax","b2d_er","b2d_ez","b2d"))
	{
		logfile<<"Error while summing "<<"b2d"<<endl;
		exit(1);
	}

	nPntE0=0;
	inf.open("xyzVectorE0");
	if (inf)
	{
		inf>>nPntE0;
		inf.close();
	}
	inf.clear();

	if(nPntE0)
	{
		vector<PointXYZ> TgCompE0;

		inf.open("TgCompE0");
		if (!inf)
		{
			logfile << "Error in open file " << "TgCompE0" << endl;
			cout << "Error in open file " << "TgCompE0" << endl;
			return 1;
		}
		TgCompE0.resize(nPntE0);
		for (i = 0; i < nPntE0; i++)
		{
			inf >> TgCompE0[i].x >> TgCompE0[i].y >> TgCompE0[i].z;
		}
		inf.close();
		inf.clear();

		if(AddDoubleFile("e_s_ax.dat","e_s_er.dat","e_s_ez.dat","e_s.dat", nPntE0,TgCompE0,npls))
		{
			logfile<<"Error while summing "<<"e_s.dat"<<endl;
			logfile.close();
			logfile.clear();
			exit(1);
		}

		if(AddDoubleFile("e_c_ax.dat","e_c_er.dat","e_c_ez.dat","e_c.dat", nPntE0,TgCompE0,npls))
		{
			logfile<<"Error while summing "<<"e_s.dat"<<endl;
			logfile.close();
			logfile.clear();
			exit(1);
		}
	}
	else
	{
		ofp.open("e_s.dat");
		ofp.close();
		ofp.clear();

		ofp.open("e_c.dat");
		ofp.close();
		ofp.clear();
	}

	RecvPlsSzB.resize(npls + 1);
	RecvPlsSzE.resize(npls + 1);

	inf.open("recvsb");
	if (!inf)
	{
		logfile << "Error in open file " << "recvsb" << endl;
		cout << "Error in open file " << "recvsb" << endl;
		return 1;
	}
	for (i = 0; i < npls; i++){ inf >> RecvPlsSzB[i];}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if (!inf)
	{
		logfile << "Error in open file " << "recvse" << endl;
		cout << "Error in open file " << "recvse" << endl;
		return 1;
	}
	for (i = 0; i < npls; i++){ inf >> RecvPlsSzE[i]; }
	inf.close();
	inf.clear();

	inf.open("b2d");
	if (!inf)
	{
		logfile << "Error in open file " << "b2d" << endl;
		cout << "Error in open file " << "b2d" << endl;
		return 1;
	}
	for (ipls = 0; ipls < npls; ipls++)
	{
		sprintf(buf, "b2d.%d", ipls + 1);
		ofp.open(buf);
		ofp << scientific << setprecision(14);
		for (i = 0; i < RecvPlsSzB[ipls]; i++)
		{
			inf >> FldHrm[0] >> FldHrm[1] >> FldHrm[2] >> FldHrm[3] >> FldHrm[4] >> FldHrm[5];
			ofp << FldHrm[0] << ' ' << FldHrm[1] << ' ' << FldHrm[2] << ' ' << FldHrm[3] << ' ' << FldHrm[4] << ' ' << FldHrm[5] << '\n';
		}
		ofp.close();
		ofp.clear();
	}
	inf.close();
	inf.clear();

	inf.open("e2d");
	if (!inf)
	{
		logfile << "Error in open file " << "e2d" << endl;
		cout << "Error in open file " << "e2d" << endl;
		return 1;
	}
	for (ipls = 0; ipls < npls; ipls++)
	{
		sprintf(buf, "e2d.%d", ipls + 1);
		ofp.open(buf);
		ofp << scientific << setprecision(14);
		for (i = 0; i < RecvPlsSzE[ipls]; i++)
		{
			inf >> FldHrm[0] >> FldHrm[1] >> FldHrm[2] >> FldHrm[3] >> FldHrm[4] >> FldHrm[5];
			ofp << FldHrm[0] << ' ' << FldHrm[1] << ' ' << FldHrm[2] << ' ' << FldHrm[3] << ' ' << FldHrm[4] << ' ' << FldHrm[5] << '\n';
		}
		ofp.close();
		ofp.clear();
	}
	inf.close();
	inf.clear();

	logfile.close();
	logfile.clear();

	return 0;
}
