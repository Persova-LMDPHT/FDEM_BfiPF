/*
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *                       FDEM_BfiPF  
 * The file contains code for run programm modules and some utility functions
 *
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024
*/

#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <direct.h>
#include <stdlib.h>
#include <string.h>
#include <windows.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

ofstream logfile;

bool isFileExist(char *file_name)
{
	ifstream inf;
	bool flag;
	flag = false;
	inf.open(file_name);
	if (inf)
	{
		flag = true;
		inf.close();
	}
	inf.clear();
	return flag;
}

int CreateProcessForEXE(char *cmdline, char *workdir)
{
	int retp;
	STARTUPINFOA si;
	PROCESS_INFORMATION pi;
	char str[256];
	ZeroMemory(&si, sizeof(si));
	si.cb = sizeof(si);
	ZeroMemory(&pi, sizeof(pi));
	if (!(retp=CreateProcessA(NULL,(LPSTR)(const char*)cmdline,NULL,NULL,FALSE,NULL,NULL,workdir,&si,&pi)))
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

void CopyFileToDir(char *CmdStr,char *fname,char *PathTo2D)
{
	cout<<"Coping file: "<<fname<<" ";
	sprintf(CmdStr,"copy %s %s\\%s",fname,PathTo2D,fname);
	system(CmdStr);
}

int main(int argc, char **argv)
{
	ifstream inf;

	logfile.open("LogNuclearHarm");

	try
	{
		ifstream inf;
		ofstream ofp;
		int i, retp, f1D, npntE0;
		char PathTo2D[256], CmdStr[256], CurDir[256], buf[1024];

		cout << "start" << endl;

		f1D = isFileExist("f2d");

		if (!f1D)
		{
			retp = CreateProcessForEXE("..\\Modules\\RegularMeshBuilder.exe", NULL);
			if (retp)
			{
				logfile << "Error: " << "RegularMeshBuilder.exe " << "returned " << retp << endl;
				logfile.close();
				exit(retp);
			}

			retp = CreateProcessForEXE("..\\Modules\\UnloadAnomalHarm.exe", NULL);
			if (retp)
			{
				logfile << "Error: " << "UnloadAnomalHarm.exe " << "returned " << retp << endl;
				logfile.close();
				exit(retp);
			}
		}
		else
		{
			ofp.open("xyzVectorE0");
			ofp << 0 << '\n';
			ofp.close();
			ofp.clear();

			ofp.open("TgCompE0");
			ofp.close();
			ofp.clear();

			ofp.open("mtr3d2d");
			ofp.close();
			ofp.clear();
		}
		_mkdir("2d");

		sprintf(PathTo2D, "2d");

		CopyFileToDir(CmdStr, "met2d", PathTo2D);
		CopyFileToDir(CmdStr, "RegularMeshBuilderSettings2D.cfg", PathTo2D);

		CopyFileToDir(CmdStr, "nthreads.txt", PathTo2D);
		CopyFileToDir(CmdStr, "numberofabdipoles", PathTo2D);
		CopyFileToDir(CmdStr, "numberofmnpoints", PathTo2D);
		CopyFileToDir(CmdStr, "xyzmn", PathTo2D);
		CopyFileToDir(CmdStr, "xyzVectorB", PathTo2D);
		CopyFileToDir(CmdStr, "xyzVectorE", PathTo2D);
		CopyFileToDir(CmdStr, "xyzVectorE0", PathTo2D);
		CopyFileToDir(CmdStr, "xyzVectorE0n", PathTo2D);
		CopyFileToDir(CmdStr, "TgCompE0", PathTo2D);
		CopyFileToDir(CmdStr, "z_sig_2d", PathTo2D);
		CopyFileToDir(CmdStr, "nu", PathTo2D);
		CopyFileToDir(CmdStr, "mlayers", PathTo2D);
		CopyFileToDir(CmdStr, "mlayersAdditional", PathTo2D);
		CopyFileToDir(CmdStr, "group", PathTo2D);
		CopyFileToDir(CmdStr, "sours", PathTo2D);
		CopyFileToDir(CmdStr, "recvsb", PathTo2D);
		CopyFileToDir(CmdStr, "recvse", PathTo2D);
		CopyFileToDir(CmdStr, "lin", PathTo2D);
		CopyFileToDir(CmdStr, "mtr3D2D", PathTo2D);
		CopyFileToDir(CmdStr, "ifreq", PathTo2D);

		CopyFileToDir(CmdStr, "marine_settingsTD", PathTo2D);
		CopyFileToDir(CmdStr, "marine_settingsFDC", PathTo2D);
		CopyFileToDir(CmdStr, "marine_settingsFD", PathTo2D);

		CopyFileToDir(CmdStr, "geoprep.dat", PathTo2D);
		CopyFileToDir(CmdStr, "3dmeshregular", PathTo2D);

		_chdir(PathTo2D);

		retp = CreateProcessForEXE("..\\..\\Modules\\CalcHarm2DHEL.exe", NULL);
		if (retp)
		{
			logfile << "Error: " << "CalcHarm2D " << "returned " << retp << endl;
			logfile.close();
			exit(retp);
		}

		system("copy e2d* ..");
		system("copy b2d* ..");
		system("copy e_s.dat ..");
		system("copy e_c.dat ..");

		_chdir("..");

		npntE0 = 0;
		inf.open("xyzVectorE0");
		if (inf)
		{
			inf >> npntE0;
			inf.close();
		}
		inf.clear();

		if (!f1D && npntE0)
		{
			ofp.open("Harm3DParams");
			ofp << 0 << '\n';
			ofp << 2 << '\n';
			ofp.close();
			ofp.clear();

			retp = CreateProcessForEXE("..\\Modules\\CalcHarm3D.exe", NULL);
			if (retp)
			{
				logfile << "Error: " << "CalcHarm3D.exe " << "returned " << retp << endl;
				logfile.close();
				exit(retp);
			}
		}

		retp = CreateProcessForEXE("..\\Modules\\SumHarm2D3D.exe", NULL);
		if (retp)
		{
			logfile << "Error: " << "SumHarm2D3D.exe " << "returned " << retp << endl;
			logfile.close();
			exit(retp);
		}

		system("del /q e2d e2d.ax e2d.u e2d_er e2d_ez e2d_ax b2d b2d.ax b2d.u b2d_er b2d_ez b2d_ax");

		system("copy e3d.* ..\\Results");
		system("copy b3d.* ..\\Results");
		system("copy e2d.* ..\\Results");
		system("copy b2d.* ..\\Results");

		system("del /s /q *.*");
		system("rmdir /s /q 2d");

		cout << "finish" << endl;
	}
	catch (...)
	{
		logfile << "Error detected" << endl;
		exit(1);
	}

	logfile.close();
	logfile.clear();

	return 0;
}
