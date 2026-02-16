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
 *  This file contains code of mesh functions
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#include "stdafx.h"
#include "task2d.h"

task2d::task2d(char *_path)
{
	strcpy(path,_path);

	rect=NULL;
	pnt=NULL;

	nmat3d=0;
	mtr3d2d=NULL;

	nmat2d=0;
	sigma=NULL;
	sigmaZ=NULL;

	reg=NULL;
	rm=NULL;
	zm=NULL;
}

task2d::~task2d()
{
	if(rect){delete [] rect; rect=NULL;}
	if(pnt){delete [] pnt; pnt=NULL;}
	if(mtr3d2d){delete [] mtr3d2d; mtr3d2d=NULL;}
	if(sigma){delete [] sigma; sigma=NULL;}
	if(sigmaZ){delete [] sigmaZ; sigmaZ=NULL;}
	if(reg){delete [] reg; reg=NULL;}
	if(rm){delete [] rm; rm=NULL;}
	if(zm){delete [] zm; zm=NULL;}
}

int task2d::Read()
{
	int i,j,k,retc;
	char buf[256];
	ifstream inf, rzf, nvtrf, nvkatf, l1f, nvk1f;
	double sum;

	sprintf(buf,"%s\\inf2tr.dat",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	inf.ignore(1000, '\n');
	inf.ignore(1000, '='); inf>>kpnt;
	inf.ignore(1000, '='); inf>>krect;
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\rz.dat",path);
	rzf.open(buf, ios::binary);
	if (!rzf) return RETCODE_NOFILE;
	pnt=new PointRZ[kpnt];
	if (!pnt) return RETCODE_NOMEM;
	for (i=0; i<kpnt; i++)
	{
		rzf > pnt[i].r;
		rzf > pnt[i].z;
		if (fabs(pnt[i].z) < 1e-10) { pnt[i].z = 1e-10; }
	}
	rzf.close();
	rzf.clear();

	sprintf(buf,"%s\\nvtr.dat",path);
	nvtrf.open(buf, ios::binary);
	if (!nvtrf) return RETCODE_NOFILE;
	rect=new Rect[krect];
	if (!rect) return RETCODE_NOMEM;
	for (i=0; i<krect; i++)
	{
		nvtrf > rect[i].nodes[2] > rect[i].nodes[3] > rect[i].nodes[0] > rect[i].nodes[1] > rect[i].nodes[4] > rect[i].rtype;
	}
	nvtrf.close();
	nvtrf.clear();

	sprintf(buf,"%s\\nvkat2d.dat",path);
	nvkatf.open(buf, ios::binary);
	if (!nvkatf) return RETCODE_NOFILE;
	for (i=0; i<krect; i++)
	{
		nvkatf > rect[i].mtr;
	}
	nvkatf.close();
	nvkatf.clear();

	nc=kpnt;

	sprintf(buf,"%s\\sigma",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	nmat2d=0;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		if(k>nmat2d)nmat2d=k;
	}
	inf.close();
	inf.clear();

	sigma=new double[nmat2d];
	sigmaZ=new double[nmat2d];
	dpr=new double[nmat2d];

	sprintf(buf,"%s\\sigma",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		sigma[k-1]=sum;
		if(sigma[k-1]<1e-12){sigma[k-1]=1e-12;}
	}
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\sigmaZ",path);
	inf.open(buf);
	if(inf)
	{
		while(!inf.eof())
		{
			inf>>k;
			if(inf.eof() || !inf.good())break;
			inf>>sum;
			sigmaZ[k-1]=sum;
			if(sigmaZ[k-1]<1e-12){sigmaZ[k-1]=1e-12;}
		}
		inf.close();
		inf.clear();
	}
	else
	{
		for(i=0;i<nmat2d;i++){sigmaZ[i]=sigma[i];}
	}

	sprintf(buf,"%s\\dpr",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>sum;
		dpr[k-1]=sum*8.84194128288307421e-12;
	}
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\mtr3d2d",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	nmat3d=0;
	while(!inf.eof())
	{
		inf>>k;
		if(inf.eof() || !inf.good())break;
		inf>>j;
		if(k>nmat3d)nmat3d=k;
	}
	inf.close();
	inf.clear();

	if(nmat3d)
	{
		mtr3d2d=new int[nmat3d];
		if (!mtr3d2d) return RETCODE_NOMEM;

		sprintf(buf,"%s\\mtr3d2d",path);
		inf.open(buf);
		if (!inf) return RETCODE_NOFILE;
		nmat3d=0;
		while(!inf.eof())
		{
			inf>>k;
			if(inf.eof() || !inf.good())break;
			inf>>j;
			mtr3d2d[k-1]=j;
		}
		inf.close();
		inf.clear();
	}

	sprintf(buf,"%s\\rz.txt",path);
	inf.open(buf);
	if (!inf) return RETCODE_NOFILE;
	inf>>qr>>qz;
	inf.close();
	inf.clear();

	nreg=(qr-1)*(qz-1);

	if(!(reg=new int[nreg])) return RETCODE_NOMEM;
	if(!(rm=new double[qr])) return RETCODE_NOMEM;
	if(!(zm=new double[qz])) return RETCODE_NOMEM;

	sprintf(buf,"%s\\r.dat",path);
	inf.open(buf, ios::binary);
	if (!inf) return RETCODE_NOFILE;
	for(i=0;i<qr;i++){inf>rm[i];}
	inf.close();
	inf.clear();

	sprintf(buf,"%s\\z.dat",path);
	inf.open(buf, ios::binary);
	if (!inf) return RETCODE_NOFILE;
	for(i=0;i<qz;i++)
	{
		inf>zm[i];
		if (fabs(zm[i]) < 1e-10) { zm[i] = 1e-10; }
	}
	inf.close();
	inf.clear();

	for(i=0;i<nreg;i++){reg[i]=i+1;}

	return RETCODE_OK;
}
