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
 *  This file contains the code for output E and B fields in receivers.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#include "stdafx.h"
#include "SmoothOutput2D.h"

bool UseRotU;
bool UseGradDivA;

struct PointXYZ
{
	double x,y,z;
};

ofstream logfile;

void AddReciver(PointXYZ &xyzVectorPi,PointXYZ &GenMin,PointRZ &recPi,double &sinfiPi,double &cosfiPi,double rc0)
{
	double x,y,r;
	x=xyzVectorPi.x-GenMin.x;
	y=xyzVectorPi.y-GenMin.y;
	r=sqrt(x*x+y*y);
	if(r<rc0)r=rc0;
	recPi.r=r;
	recPi.z=xyzVectorPi.z;
	sinfiPi=y/r;
	cosfiPi=x/r;
}

struct Line
{
	PointXYZ A,B;
};

void Unpuk(double *uq,int npls,int slae_n,int tmap_n)
{
	int i,ipls;
	for(ipls=npls-1;ipls>0;ipls--)
	{
		for(i=0;i<slae_n;i++)
		{
			uq[ipls*tmap_n+slae_n-i-1]=uq[(ipls+1)*slae_n-i-1];
		}
	}
}

wchar_t *convertCharArrayToLPCWSTR(const char* charArray)
{
	wchar_t* wString = new wchar_t[4096];
	MultiByteToWideChar(CP_ACP, 0, charArray, -1, wString, 4096);
	return wString;
}

int main()
{
	int i,k,l,retp,nthreads;
	ifstream inf;
	ofstream ofp,ofps,ofpc;
	double *U[2],*B[2],*Er[2],*Ez[2],*Er0[2],*Ez0[2],*Ax[2],nu,w;
	FILE *fp;
	double x,y,z;
	SmoothOutput2D SobjU,SobjB,SobjE0;

	PointXYZ *xyzVectorE,*xyzVectorB,*xyzVectorE0;
	PointXYZ Av,Bv,Tp;

	int npntE,npntB,npntE0,npntE0or;
	int *matrec,*matrec0;
	PointRZ *recE,*recB,*recE0;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double *sinfiE0,*cosfiE0;
	double DminSrc;

	int p1,p2,ipls,npls,nprof;
	vector<Line> GenLin;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecvPlsIgE0, RecToSourceB, RecToSourceE, RecToSourceE0;

	logfile.open("LogOutput2d_AVU");

	nthreads = 1;
	inf.open("nthreads.txt");
	if (inf)
	{
		inf >> nthreads;
		inf.close();
	}
	inf.clear();

	if (nthreads < 1) { nthreads = 1; }

	omp_set_num_threads(nthreads);

	npls=0;
	inf.open("group");
	if(!inf)
	{
		logfile<<"Error in open file "<<"group"<<endl;
		cout<<"Error in open file "<<"group"<<endl;
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

	RecvPlsIgB.resize(npls+1);
	RecvPlsIgE.resize(npls+1);
	RecvPlsIgE0.resize(npls+1);

	inf.open("xyzVectorE0");
	if(inf)
	{
		inf>>npntE0;
		inf.close();
	}
	else
	{
		npntE0=0;
	}
	inf.clear();

	inf.open("recvsb");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvsb"<<endl;
		cout<<"Error in open file "<<"recvsb"<<endl;
		return 1;
	}
	RecvPlsIgB[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgB[i+1];}
	inf.close();
	inf.clear();

	inf.open("recvse");
	if(!inf)
	{
		logfile<<"Error in open file "<<"recvse"<<endl;
		cout<<"Error in open file "<<"recvse"<<endl;
		return 1;
	}
	RecvPlsIgE[0]=0;
	for(i=0;i<npls;i++){inf>>RecvPlsIgE[i+1];}
	inf.close();
	inf.clear();

	RecvPlsIgE0[0]=0;
	for(i=0;i<npls;i++){RecvPlsIgE0[i+1]=npntE0;}

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
		RecvPlsIgE0[i+1]+=RecvPlsIgE0[i];
	}

	npntB=RecvPlsIgB[npls];
	npntE=RecvPlsIgE[npls];
	npntE0=RecvPlsIgE0[npls];

	RecToSourceB.resize(npntB);
	RecToSourceE.resize(npntE);
	RecToSourceE0.resize(npntE0);

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgB[i];k<RecvPlsIgB[i+1];k++)
		{
			RecToSourceB[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE[i];k<RecvPlsIgE[i+1];k++)
		{
			RecToSourceE[k]=i;
		}
	}

	for(i=0;i<npls;i++)
	{
		for(k=RecvPlsIgE0[i];k<RecvPlsIgE0[i+1];k++)
		{
			RecToSourceE0[k]=i;
		}
	}

	GenLin.resize(npls);
	inf.open("sours");
	if(!inf)
	{
		logfile<<"Error in open file "<<"sours"<<endl;
		cout<<"Error in open file "<<"sours"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenLin[i].A.x>>GenLin[i].A.y>>GenLin[i].A.z>>GenLin[i].B.x>>GenLin[i].B.y>>GenLin[i].B.z;}
	inf.close();
	inf.clear();

	task2d task("Ax");
	double rc0;

	retp=task.Read();
	if(retp)
	{
		logfile<<"Function Read() for Ax returned "<<retp<<'\n';
		return 1;
	}

	rc0=0.5*(task.rm[0]+task.rm[1]);

	DminSrc=rc0;
	inf.open("dminsrc");
	if(inf)
	{
		inf>>DminSrc;
		inf.close();
	}
	inf.clear();

	inf.open("xyzVectorB");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorB"<<'\n';
		return 1;
	}
	inf>>npntB;
	xyzVectorB=new PointXYZ[npntB];
	recB=new PointRZ[npntB*2];
	if(!recB)return RETCODE_NOMEM;
	sinfiB=new double[npntB*2];
	if(!sinfiB)return RETCODE_NOMEM;
	cosfiB=new double[npntB*2];
	if(!cosfiB)return RETCODE_NOMEM;
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
		{
			Av=GenLin[RecToSourceB[i]].A;
			Bv=GenLin[RecToSourceB[i]].B;
			inf>>xyzVectorB[i].x>>xyzVectorB[i].y>>xyzVectorB[i].z;
			AddReciver(xyzVectorB[i],Av,recB[2*i],sinfiB[2*i],cosfiB[2*i],rc0);
			AddReciver(xyzVectorB[i],Bv,recB[2*i+1],sinfiB[2*i+1],cosfiB[2*i+1],rc0);
		}
	}
	inf.close();
	inf.clear();

	inf.open("xyzVectorE");
	if(!inf)
	{
		logfile<<"No file "<<"xyzVectorE"<<'\n';
		return 1;
	}
	inf>>npntE;
	xyzVectorE=new PointXYZ[npntE];
	recE=new PointRZ[npntE*2];
	if(!recE)return RETCODE_NOMEM;
	sinfiE=new double[npntE*2];
	if(!sinfiE)return RETCODE_NOMEM;
	cosfiE=new double[npntE*2];
	if(!cosfiE)return RETCODE_NOMEM;
	matrec=new int[npntE*2];
	if(!matrec)return RETCODE_NOMEM;
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=RecvPlsIgE[ipls];i<RecvPlsIgE[ipls+1];i++)
		{
			Av=GenLin[RecToSourceE[i]].A;
			Bv=GenLin[RecToSourceE[i]].B;
			inf>>xyzVectorE[i].x>>xyzVectorE[i].y>>xyzVectorE[i].z;
			AddReciver(xyzVectorE[i],Av,recE[2*i],sinfiE[2*i],cosfiE[2*i],rc0);
			AddReciver(xyzVectorE[i],Bv,recE[2*i+1],sinfiE[2*i+1],cosfiE[2*i+1],rc0);
			matrec[2*i]=0;
			matrec[2*i+1]=0;
		}
	}
	inf.close();
	inf.clear();

	npntE0or=0;
	if(npntE0)
	{
		inf.open("xyzVectorE0");
		if(!inf)
		{
			logfile<<"No file "<<"xyzVectorE0"<<'\n';
			return 1;
		}
		inf>>npntE0or;
		xyzVectorE0=new PointXYZ[npntE0or];
		for(i=0;i<npntE0or;i++)
		{
			inf>>Tp.x>>Tp.y>>Tp.z;
			xyzVectorE0[i]=Tp;
		}
		inf.close();
		inf.clear();
	}

	recE0=new PointRZ[npntE0*2];
	if(!recE0)return RETCODE_NOMEM;
	sinfiE0=new double[npntE0*2];
	if(!sinfiE0)return RETCODE_NOMEM;
	cosfiE0=new double[npntE0*2];
	if(!cosfiE0)return RETCODE_NOMEM;
	matrec0=new int[npntE0*2];
	if(!matrec0)return RETCODE_NOMEM;
	for(ipls=0;ipls<npls;ipls++)
	{
		l=0;
		for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
		{
			Av=GenLin[RecToSourceE0[i]].A;
			Bv=GenLin[RecToSourceE0[i]].B;
			Tp=xyzVectorE0[l];
			l++;
			AddReciver(Tp,Av,recE0[2*i],sinfiE0[2*i],cosfiE0[2*i],rc0);
			AddReciver(Tp,Bv,recE0[2*i+1],sinfiE0[2*i+1],cosfiE0[2*i+1],rc0);
			matrec0[2*i]=0;
			matrec0[2*i+1]=0;
		}
	}

	if(npntE0)
	{
		inf.open("xyzVectorE0n");
		if(inf)
		{
			for(i=0;i<npntE0or;i++)
			{
				inf>>matrec0[2*i];
				matrec0[2*i+1]=matrec0[2*i];
			}
			for(ipls=1;ipls<npls;ipls++)
			{
				for(i=0;i<npntE0or;i++)
				{
					matrec0[2*(RecvPlsIgE0[ipls]+i)]=matrec0[2*i];
					matrec0[2*(RecvPlsIgE0[ipls]+i)+1]=matrec0[2*i+1];
				}
			}
			inf.close();

			for(ipls=0;ipls<npls;ipls++)
			{
				for(i=RecvPlsIgE0[ipls+1];i<RecvPlsIgE0[ipls+1];i++)
				{
					matrec0[2*i]=task.mtr3d2d[matrec0[2*i]-1];
					matrec0[2*i+1]=matrec0[2*i];
				}
			}
		}
		inf.clear();
	}


	inf.open("nu");
	if(!inf)
	{
		logfile<<"No file "<<"nu"<<'\n';
		return 1;
	}
	inf>>nu;
	inf.close();
	inf.clear();

	w=2.0*PI*nu;

	U[0]=new double[task.kpnt*npls];
	if(!U[0])return RETCODE_NOMEM;
	fp=fopen("u.re","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(U[0],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	U[1]=new double[task.kpnt*npls];
	if(!U[1])return RETCODE_NOMEM;
	fp=fopen("u.im","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(U[1],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	Ax[0]=new double[task.kpnt*npls];
	if(!Ax[0])return RETCODE_NOMEM;
	fp=fopen("Ax\\v2s.dat","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(Ax[0],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	Ax[1]=new double[task.kpnt*npls];
	if(!Ax[1])return RETCODE_NOMEM;
	fp=fopen("Ax\\v2c.dat","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(Ax[1],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	int MemSize;
	HANDLE hMapFileSin, hMapFileCos;
	LPCTSTR pBufSin, pBufCos;
	double *BuffSin, *BuffCos;
	int irecs, irecc;
	int ifreq;
	char ESinName[256], ECosName[256];

	BuffSin=BuffCos=NULL;

	if(npntE0or)
	{
		ifreq=0;
		inf.open("ifreq");
		if (!inf)
		{
			cout << "Error in open file " << "ifreq" << endl;
			logfile << "Error in open file " << "ifreq" << endl;
			return 1;
		}
		inf>>ifreq;
		inf.close();
		inf.clear();

		sprintf(ESinName,"Esin_%d",ifreq);
		sprintf(ECosName,"Ecos_%d",ifreq);

		MemSize = npntE0or * 3 * npls * sizeof(double);

		BuffSin = new double[npntE0or * 3 * npls];
		BuffCos = new double[npntE0or * 3 * npls];

		cout << "Memmory for E norm " << (2.0 * MemSize) / (1024 * 1024) << "Mb" << '\n';

		hMapFileSin = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(ESinName));
		if (hMapFileSin == NULL || hMapFileSin == INVALID_HANDLE_VALUE)
		{
			logfile<<"Cannot create a object in memory ("<<GetLastError()<<")."<<endl;
			return 1;
		}
		pBufSin = (LPTSTR)MapViewOfFile(hMapFileSin, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBufSin == NULL)
		{
			logfile<<"File in memory cannot be represented ("<<GetLastError()<<")."<<endl;
			return 1;
		}

		hMapFileCos = CreateFileMapping(INVALID_HANDLE_VALUE, NULL, PAGE_READWRITE, 0, MemSize, convertCharArrayToLPCWSTR(ECosName));
		if (hMapFileCos == NULL || hMapFileCos == INVALID_HANDLE_VALUE)
		{
			logfile<<"Cannot create a object in memory ("<<GetLastError()<<")."<<endl;
			return 1;
		}
		pBufCos = (LPTSTR)MapViewOfFile(hMapFileCos, FILE_MAP_ALL_ACCESS, 0, 0, MemSize);
		if (pBufCos == NULL)
		{
			logfile<<"File in memory cannot be represented ("<<GetLastError()<<")."<<endl;
			return 1;
		}

		CopyMemory((char *)BuffSin, (PVOID)pBufSin, MemSize);
		CopyMemory((char *)BuffCos, (PVOID)pBufCos, MemSize);
	}


	retp=SobjU.Init(npntE*2,task.kpnt,task.krect,recE,task.pnt,task.rect,0,NULL,task.sigma,task.sigmaZ,task.dpr,task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm,w);
	if(retp)
	{
		logfile<<"SobjU.Init() returned "<<retp<<endl;
		return retp;
	}

	retp=SobjB.Init(npntB*2,task.kpnt,task.krect,recB,task.pnt,task.rect,0,NULL,task.sigma,task.sigmaZ,task.dpr,task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm,w);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp<<endl;
		return retp;
	}

	retp=SobjE0.Init(npntE0*2,task.kpnt,task.krect,recE0,task.pnt,task.rect,0,NULL,task.sigma,task.sigmaZ,task.dpr,task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm,w);
	if(retp)
	{
		logfile<<"SobjE0.Init() returned "<<retp<<endl;
		return retp;
	}

	B[0]=new double[npntB*2];
	if(!B[0])return RETCODE_NOMEM;
	B[1]=new double[npntB*2];
	if(!B[1])return RETCODE_NOMEM;

	Er[0]=new double[npntE*2];
	if(!Er[0])return RETCODE_NOMEM;
	Er[1]=new double[npntE*2];
	if(!Er[1])return RETCODE_NOMEM;

	Ez[0]=new double[npntE*2];
	if(!Ez[0])return RETCODE_NOMEM;
	Ez[1]=new double[npntE*2];
	if(!Ez[1])return RETCODE_NOMEM;

	Er0[0]=new double[npntE0*2];
	if(!Er0[0])return RETCODE_NOMEM;
	Er0[1]=new double[npntE0*2];
	if(!Er0[1])return RETCODE_NOMEM;

	Ez0[0]=new double[npntE0*2];
	if(!Ez0[0])return RETCODE_NOMEM;
	Ez0[1]=new double[npntE0*2];
	if(!Ez0[1])return RETCODE_NOMEM;

	k=npntE*2;
	for(i=0;i<k;i++){Er[0][i]=Er[1][i]=Ez[0][i]=Ez[1][i]=0.0;}

	k=npntE0*2;
	for(i=0;i<k;i++){Er0[0][i]=Er0[1][i]=Ez0[0][i]=Ez0[1][i]=0.0;}

	UseRotU=true;
	UseGradDivA=true;


	SobjU.Output_U_Er(U[0],U[1],Ax[0],Ax[1],Er[0],Er[1],1,RecToSourceE);
	SobjU.Output_U_Ez(U[0],U[1],Ax[0],Ax[1],Ez[0],Ez[1],1,RecToSourceE);

	SobjE0.Output_U_Er(U[0],U[1],Ax[0],Ax[1],Er0[0],Er0[1],1,RecToSourceE0);
	SobjE0.Output_U_Ez(U[0],U[1],Ax[0],Ax[1],Ez0[0],Ez0[1],1,RecToSourceE0);

	SobjB.Output_U_B(U[0],U[1],B[0],B[1],1,RecToSourceB);


	ofp.open("e2d.u");
	ofp<<scientific<<setprecision(14);
	for(i=0;i<npntE;i++)
	{
		x=(Er[0][2*i+1]*cosfiE[2*i+1]-Er[0][2*i]*cosfiE[2*i]);
		y=(Er[0][2*i+1]*sinfiE[2*i+1]-Er[0][2*i]*sinfiE[2*i]);
		z=(Ez[0][2*i+1]-Ez[0][2*i]);
		ofp<<-x<<' '<<-y<<' '<<-z<<' ';
		x=(Er[1][2*i+1]*cosfiE[2*i+1]-Er[1][2*i]*cosfiE[2*i]);
		y=(Er[1][2*i+1]*sinfiE[2*i+1]-Er[1][2*i]*sinfiE[2*i]);
		z=(Ez[1][2*i+1]-Ez[1][2*i]);
		ofp<<-x<<' '<<-y<<' '<<-z<<'\n';
	}
	ofp.close();
	ofp.clear();


	irecs=irecc=0;
	for(ipls=0;ipls<npls;ipls++)
	{
		for(i=RecvPlsIgE0[ipls];i<RecvPlsIgE0[ipls+1];i++)
		{
			if(recE0[i].r>DminSrc || (fabs(recE0[i].z-GenLin[ipls].A.z)>DminSrc && fabs(recE0[i].z-GenLin[ipls].B.z)>DminSrc))
			{
				x=(Er0[0][2*i+1]*cosfiE0[2*i+1]-Er0[0][2*i]*cosfiE0[2*i]);
				y=(Er0[0][2*i+1]*sinfiE0[2*i+1]-Er0[0][2*i]*sinfiE0[2*i]);
				z=(Ez0[0][2*i+1]-Ez0[0][2*i]);
				BuffSin[irecs]+=-x; irecs++;
				BuffSin[irecs]+=-y; irecs++;
				BuffSin[irecs]+=-z; irecs++;
					
				x=(Er0[1][2*i+1]*cosfiE0[2*i+1]-Er0[1][2*i]*cosfiE0[2*i]);
				y=(Er0[1][2*i+1]*sinfiE0[2*i+1]-Er0[1][2*i]*sinfiE0[2*i]);
				z=(Ez0[1][2*i+1]-Ez0[1][2*i]);
				BuffCos[irecc]+=-x; irecc++;
				BuffCos[irecc]+=-y; irecc++;
				BuffCos[irecc]+=-z; irecc++;
			}
			else
			{
				BuffSin[irecs]+=0.0; irecs++;
				BuffSin[irecs]+=0.0; irecs++;
				BuffSin[irecs]+=0.0; irecs++;
					
				BuffCos[irecc]+=0.0; irecc++;
				BuffCos[irecc]+=0.0; irecc++;
				BuffCos[irecc]+=0.0; irecc++;
			}
		}
	}

	if(npntE0or)
	{
		CopyMemory((PVOID)pBufSin,(char *)BuffSin, MemSize);
		CopyMemory((PVOID)pBufCos,(char *)BuffCos, MemSize);

		UnmapViewOfFile(pBufSin);
		UnmapViewOfFile(pBufCos);

		CloseHandle(hMapFileSin);
		CloseHandle(hMapFileCos);
	}


	ofp.open("b2d.u");
	ofp<<scientific<<setprecision(14);
	for(i=0;i<npntB;i++)
	{
		x=(-B[0][2*i+1]*sinfiB[2*i+1]+B[0][2*i]*sinfiB[2*i]);
		y=(B[0][2*i+1]*cosfiB[2*i+1]-B[0][2*i]*cosfiB[2*i]);
		ofp<<-x<<' '<<-y<<' '<<0.0<<' ';
		x=(-B[1][2*i+1]*sinfiB[2*i+1]+B[1][2*i]*sinfiB[2*i]);
		y=(B[1][2*i+1]*cosfiB[2*i+1]-B[1][2*i]*cosfiB[2*i]);
		ofp<<-x<<' '<<-y<<' '<<0.0<<'\n';
	}
	ofp.close();
	ofp.clear();

	if(U[0]){delete [] U[0]; U[0]=NULL;}
	if(U[1]){delete [] U[1]; U[1]=NULL;}

	if(B[0]){delete [] B[0]; B[0]=NULL;}
	if(B[1]){delete [] B[1]; B[1]=NULL;}

	if(Ax[0]){delete [] Ax[0]; Ax[0]=NULL;}
	if(Ax[1]){delete [] Ax[1]; Ax[1]=NULL;}

	if(Er[0]){delete [] Er[0]; Er[0]=NULL;}
	if(Er[1]){delete [] Er[1]; Er[1]=NULL;}

	if(Ez[0]){delete [] Ez[0]; Ez[0]=NULL;}
	if(Ez[1]){delete [] Ez[1]; Ez[1]=NULL;}

	if(Er0[0]){delete [] Er[0]; Er[0]=NULL;}
	if(Er0[1]){delete [] Er[1]; Er[1]=NULL;}

	if(Ez0[0]){delete [] Ez[0]; Ez[0]=NULL;}
	if(Ez0[1]){delete [] Ez[1]; Ez[1]=NULL;}

	if(recE){delete [] recE; recE=NULL;}
	if(sinfiE){delete [] sinfiE; sinfiE=NULL;}
	if(cosfiE){delete [] cosfiE; cosfiE=NULL;}

	if(recE0){delete [] recE0; recE0=NULL;}
	if(sinfiE0){delete [] sinfiE0; sinfiE0=NULL;}
	if(cosfiE0){delete [] cosfiE0; cosfiE0=NULL;}

	if(BuffSin){delete [] BuffSin; BuffSin=NULL;}
	if(BuffCos){delete [] BuffCos; BuffCos=NULL;}

	logfile.close();
	logfile.clear();

	return 0;
}
