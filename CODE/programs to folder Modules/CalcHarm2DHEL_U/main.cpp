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
	void read(ifstream &inf)
	{
		inf>>x>>y>>z;
	}
	void write(ofstream &ofp)
	{
		ofp<<x<<' '<<y<<' '<<z;
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
		logfile<<"Error: "<<cmd<<" returned "<<retp<<endl;
		logfile.close();
		exit(retp);
	}
}

wchar_t* convertCharArrayToLPCWSTR(const char* charArray)
{
	wchar_t* wString = new wchar_t[4096];
	MultiByteToWideChar(CP_ACP, 0, charArray, -1, wString, 4096);
	return wString;
}

int SumHarmResFile(char* FileName0, char* FileName1, char* FileName2, double cf)
{
	ifstream inf1, inf2;
	ofstream ofp;
	int i;
	double t1, t2;
	inf1.open(FileName0);
	if (!inf1)return 1;
	inf2.open(FileName1);
	if (!inf2)return 1;
	ofp.open(FileName2);
	ofp << scientific << setprecision(14);
	i = 0;
	while (!inf1.eof() && !inf2.eof())
	{
		inf1 >> t1;
		if (inf1.eof() || !inf1.good())break;
		inf2 >> t2;
		if (inf2.eof() || !inf2.good())break;
		ofp << t1 + cf * t2;
		i++;
		if (i % 6)
		{
			ofp << ' ';
		}
		else
		{
			ofp << '\n';
		}
	}
	if (i % 6)return 1;
	inf1.close();
	inf1.clear();
	inf2.close();
	inf2.clear();
	ofp.close();
	ofp.clear();
	return 0;
}

int main(int argc, char **argv)
{
	ifstream inf;
	ofstream ofp;
	char ESinName[256], ECosName[256], str[256];
	HANDLE hMapFileSin, hMapFileCos;
	LPCTSTR pBufSin, pBufCos;
	int i,j,ipls,npls,nprof,k,p1,p2,npntE0,ifreq,MemSize;
	double v[6];

	logfile.open("logCalcHarm2DHEL_U");

	npls = 0;
	inf.open("group");
	if (!inf)
	{
		logfile << "Error in open file " << "group" << endl;
		cout << "Error in open file " << "group" << endl;
		return 1;
	}
	inf>>nprof;
	for(i=0;i<nprof;i++)
	{
		inf>>k>>p1>>p2;
		npls+=p2-p1+1;
	}
	inf.close();
	inf.clear();

	logfile<<"start Ax"<<endl;
	run("..\\..\\Modules\\CalcHarm2D_Ax.exe");
	logfile<<"finish Ax"<<endl;

	logfile<<"start HEL"<<endl;
	run("..\\..\\Modules\\bound.exe");
	run("..\\..\\Modules\\u.exe");
	logfile<<"finish HEL"<<endl;

	npntE0 = 0;
	inf.open("xyzVectorE0");
	if (!inf)
	{
		cout << "Error in open file " << "xyzVectorE0" << endl;
		logfile << "Error in open file " << "xyzVectorE0" << endl;
		return 1;
	}
	inf >> npntE0;
	inf.close();
	inf.clear();

	ifreq = 0;
	inf.open("ifreq");
	if (!inf)
	{
		cout << "Error in open file " << "ifreq" << endl;
		logfile << "Error in open file " << "ifreq" << endl;
		return 1;
	}
	inf >> ifreq;
	inf.close();
	inf.clear();

	if (npntE0)
	{
		MemSize = npntE0 * 3 * npls * sizeof(double);

		sprintf(ESinName, "Esin_%d", ifreq);
		sprintf(ECosName, "Ecos_%d", ifreq);

		cout << "Memmory for E norm " << (2.0 * MemSize) / (1024 * 1024) << "Mb" << '\n';

		hMapFileSin = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(ESinName));
		if (hMapFileSin == NULL || hMapFileSin == INVALID_HANDLE_VALUE)
		{
			logfile << "Cannot create a object in memory (" << GetLastError() << ")." << endl;
			return 1;
		}
		pBufSin = (LPTSTR)MapViewOfFile(hMapFileSin, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBufSin == NULL)
		{
			logfile << "File in memory cannot be represented (" << GetLastError() << ")." << endl;
			return 1;
		}

		hMapFileCos = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(ECosName));
		if (hMapFileCos == NULL || hMapFileCos == INVALID_HANDLE_VALUE)
		{
			logfile << "Cannot create a object in memory (" << GetLastError() << ")." << endl;
			return 1;
		}
		pBufCos = (LPTSTR)MapViewOfFile(hMapFileCos, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBufCos == NULL)
		{
			logfile << "File in memory cannot be represented (" << GetLastError() << ")." << endl;
			return 1;
		}
	}

	run("..\\..\\Modules\\OutputSmooth2DHarm_Ax.exe");
	run("..\\..\\Modules\\OutputSmooth2DHarm_U.exe");

	if (SumHarmResFile("e2d.ax", "e2d.u", "e2d", 1.0))
	{
		logfile << "Error while summing " << "e2d" << endl;
		exit(1);
	}

	if (SumHarmResFile("b2d.ax", "b2d.u", "b2d", 1.0))
	{
		logfile << "Error while summing " << "b2d" << endl;
		exit(1);
	}

	if (npntE0)
	{
		int nnn, fff;
		vector<PointXYZ> TgCompE0;
		float tmpf;
		PointXYZ E0;
		const int size_f = sizeof(float);
		const int size_d = sizeof(double);
		const int size_d2 = size_d / 2;

		inf.open("TgCompE0");
		if (!inf)
		{
			logfile << "Error in open file " << "TgCompE0" << endl;
			cout << "Error in open file " << "TgCompE0" << endl;
			return 1;
		}
		TgCompE0.resize(npntE0);
		for (i = 0; i < npntE0; i++)
		{
			inf >> TgCompE0[i].x >> TgCompE0[i].y >> TgCompE0[i].z;
		}
		inf.close();
		inf.clear();

		ofp.open("e_s.dat", ios::binary);
		nnn = 0;
		for (ipls = 0; ipls < npls; ipls++)
		{
			for (i = 0; i < npntE0; i++)
			{
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.x) + fff) = pBufSin[nnn + fff]; }
				nnn += size_d2;
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.y) + fff) = pBufSin[nnn + fff]; }
				nnn += size_d2;
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.z) + fff) = pBufSin[nnn + fff]; }
				nnn += size_d2;
				tmpf = (float)(E0.x * TgCompE0[i].x + E0.y * TgCompE0[i].y + E0.z * TgCompE0[i].z);
				ofp.write((char*)&tmpf, size_f);
			}
		}
		ofp.close();
		ofp.clear();

		ofp.open("e_c.dat", ios::binary);
		nnn = 0;
		for (ipls = 0; ipls < npls; ipls++)
		{
			for (i = 0; i < npntE0; i++)
			{
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.x) + fff) = pBufCos[nnn + fff]; }
				nnn += size_d2;
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.y) + fff) = pBufCos[nnn + fff]; }
				nnn += size_d2;
				for (fff = 0; fff < size_d2; fff++) { *(((wchar_t*)&E0.z) + fff) = pBufCos[nnn + fff]; }
				nnn += size_d2;
				tmpf = (float)(E0.x * TgCompE0[i].x + E0.y * TgCompE0[i].y + E0.z * TgCompE0[i].z);
				ofp.write((char*)&tmpf, size_f);
			}
		}
		ofp.close();
		ofp.clear();

		UnmapViewOfFile(pBufSin);
		UnmapViewOfFile(pBufCos);

		CloseHandle(hMapFileSin);
		CloseHandle(hMapFileCos);
	}
	else
	{
		ofp.open("e_s.dat", ios::binary);
		ofp.close();
		ofp.clear();

		ofp.open("e_c.dat", ios::binary);
		ofp.close();
		ofp.clear();
	}

	vector<int> RecvPlsIgB, RecvPlsIgE;
	vector<int> RecToSourceB, RecToSourceE;

	RecvPlsIgB.resize(npls + 1);
	RecvPlsIgE.resize(npls + 1);

	inf.open("recvsb");
	if (!inf)
	{
		logfile << "Error in open file " << "recvsb" << endl;
		cout << "Error in open file " << "recvsb" << endl;
		return 1;
	}
	RecvPlsIgB[0] = 0;
	for (i = 0; i < npls; i++) { inf >> RecvPlsIgB[i + 1]; }
	inf.close();
	inf.clear();

	inf.open("recvse");
	if (!inf)
	{
		logfile << "Error in open file " << "recvse" << endl;
		cout << "Error in open file " << "recvse" << endl;
		return 1;
	}
	RecvPlsIgE[0] = 0;
	for (i = 0; i < npls; i++) { inf >> RecvPlsIgE[i + 1]; }
	inf.close();
	inf.clear();

	for (i = 0; i < npls; i++)
	{
		RecvPlsIgB[i + 1] += RecvPlsIgB[i];
		RecvPlsIgE[i + 1] += RecvPlsIgE[i];
	}

	inf.open("e2d");
	if (!inf)
	{
		logfile << "Error in open file " << "e2d" << endl;
		cout << "Error in open file " << "e2d" << endl;
		return 1;
	}
	for (ipls = 0; ipls < npls; ipls++)
	{
		sprintf(str, "e2d.%d", ipls + 1);
		ofp.open(str);
		ofp << scientific << setprecision(14);
		for (i = RecvPlsIgE[ipls]; i < RecvPlsIgE[ipls + 1]; i++)
		{
			inf >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5];
			ofp << v[0] << '\t' << v[1] << '\t' << v[2] << '\t' << v[3] << '\t' << v[4] << '\t' << v[5] << '\n';
		}
		ofp.close();
		ofp.clear();
	}
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
		sprintf(str, "b2d.%d", ipls + 1);
		ofp.open(str);
		ofp << scientific << setprecision(14);
		for (i = RecvPlsIgB[ipls]; i < RecvPlsIgB[ipls + 1]; i++)
		{
			inf >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5];
			ofp << v[0] << '\t' << v[1] << '\t' << v[2] << '\t' << v[3] << '\t' << v[4] << '\t' << v[5] << '\n';
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
