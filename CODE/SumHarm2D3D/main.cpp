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
 *  This file contains code for summing 2D and 3D harmonic tasks results
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024 
*/

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <direct.h>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

ofstream logfile;

char RunPath[1024],str[1024];

int GetNumberOfPlaces(int& npls)
{
	int i, k, p1, p2, nprof;
	ifstream inf;

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

	return 0;
}

struct Field
{
	double Fs[3],Fc[3];
	void Clear()
	{
		Fs[0]=Fs[1]=Fs[2]=0.0;
		Fc[0]=Fc[1]=Fc[2]=0.0;
	}
	void AddField(Field &f,double cff)
	{
		Fs[0]+=f.Fs[0]*cff;
		Fs[1]+=f.Fs[1]*cff;
		Fs[2]+=f.Fs[2]*cff;
		Fc[0]+=f.Fc[0]*cff;
		Fc[1]+=f.Fc[1]*cff;
		Fc[2]+=f.Fc[2]*cff;
	}
	void Read(ifstream &inf)
	{
		inf>>Fs[0]>>Fs[1]>>Fs[2]>>Fc[0]>>Fc[1]>>Fc[2];
	}
	void Write(ofstream &ofp)
	{
		ofp<<Fs[0]<<' '<<Fs[1]<<' '<<Fs[2]<<' '<<Fc[0]<<' '<<Fc[1]<<' '<<Fc[2]<<'\n';
	}
	void SummValues(Field &f1,Field &f2)
	{
		Fs[0]=f1.Fs[0]+f2.Fs[0];
		Fs[1]=f1.Fs[1]+f2.Fs[1];
		Fs[2]=f1.Fs[2]+f2.Fs[2];
		Fc[0]=f1.Fc[0]+f2.Fc[0];
		Fc[1]=f1.Fc[1]+f2.Fc[1];
		Fc[2]=f1.Fc[2]+f2.Fc[2];
	}
};

int ReadFields(int nrec,vector<Field> &fv,char *fname,bool fmem)
{
	int irec;
	ifstream inf;
	inf.open(fname);
	if(!inf)
	{
		cout<<"Can't ope file "<<fname<<endl;
		logfile<<"Can't ope file "<<fname<<endl;
		return 1;
	}
	if(fmem)
	{
		fv.resize(nrec);
	}
	for(irec=0;irec<nrec;irec++)
	{
		fv[irec].Read(inf);
	}
	inf.close();
	inf.clear();
	return 0;
}

void WriteFields(int nrec,vector<Field> &fv,char *fname)
{
	int irec;
	ofstream ofp;
	ofp.open(fname);
	for(irec=0;irec<nrec;irec++)
	{
		fv[irec].Write(ofp);
	}
	ofp.close();
	ofp.clear();
}

void WriteFields(int beg,int end,vector<Field> &fv,char *fname)
{
	int irec;
	ofstream ofp;
	ofp.open(fname);
	ofp<<scientific<<setprecision(14);
	for(irec=beg;irec<end;irec++)
	{
		fv[irec].Write(ofp);
	}
	ofp.close();
	ofp.clear();
}

int main(int argc,char **argv)
{
	int retp,i;

	ifstream inf;
	ofstream ofp;

	int ipls,npls,nall,npntE0;

	int nRec;
	vector<Field> fnrm,fanm,fsum;
	vector<int> vNrec;

	logfile.open("SumHarm2D3D.log");

	_getcwd(RunPath, 1023);
	cout << "Starting programm in " << RunPath << '\n';
	logfile << "Starting programm in " << RunPath << '\n';

	retp=GetNumberOfPlaces(npls);
	if(retp)
	{
		cout << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		logfile << "Function GetNumberOfPlaces retrurned " << retp << '\n';
		return 1;
	}

	npntE0=0;
	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	inf.clear();

	vNrec.resize(npls);

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorB"<<endl;
		cout<<"Error in open file "<<"xyzVectorB"<<endl;
		return 1;
	}
	inf>>nRec;
	inf.close();
	inf.clear();

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		if(npntE0)
		{
			nall=vNrec[ipls];
			sprintf(str,"b2d.%d",ipls+1);
			ReadFields(nall,fnrm,str,true);
			sprintf(str,"b3d_anom.%d",ipls+1);
			ReadFields(nall,fanm,str,true);
			fsum.clear();
			fsum.resize(nall);
			for(i=0;i<nall;i++)
			{
				fsum[i].Clear();
				fsum[i].SummValues(fnrm[i],fanm[i]);
			}
			sprintf(str,"b3d.%d",ipls+1);
			WriteFields(0,nall,fsum,str);
		}
		else
		{
			sprintf(str,"copy b2d.%d b3d.%d",ipls+1,ipls+1);
			system(str);
		}
	}

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"Error in open file "<<"xyzVectorE"<<endl;
		cout<<"Error in open file "<<"xyzVectorE"<<endl;
		return 1;
	}
	inf>>nRec;
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>vNrec[i];}
	inf.close();
	inf.clear();

	for(ipls=0;ipls<npls;ipls++)
	{
		if(npntE0)
		{
			nall=vNrec[ipls];
			sprintf(str,"e2d.%d",ipls+1);
			ReadFields(nall,fnrm,str,true);
			sprintf(str,"e3d_anom.%d",ipls+1);
			ReadFields(nall,fanm,str,true);
			fsum.clear();
			fsum.resize(nall);
			for(i=0;i<nall;i++)
			{
				fsum[i].Clear();
				fsum[i].SummValues(fnrm[i],fanm[i]);
			}
			sprintf(str,"e3d.%d",ipls+1);
			WriteFields(0,nall,fsum,str);
		}
		else
		{
			sprintf(str,"copy e2d.%d e3d.%d",ipls+1,ipls+1);
			system(str);
		}
	}

	logfile.close();
	logfile.clear();

	return 0;
}
