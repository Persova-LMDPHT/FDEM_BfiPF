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
 *  This file contains code for output E and B fields in receivers.
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
#include "SimpleOutput2D.h"

struct PointXYZ
{
	double x,y,z;
};

struct SqLoop
{
	PointXYZ A,B;
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

void get_int(char *str, int &ipos, int &var)
{
	int l;
	bool flag;
	char ch;
	char word[32];
	flag = true;
	l = 0;
	do
	{
		ch = str[ipos];
		if (!(ch == ' ' || ch == '\t' || ch == '\n'))
		{
			word[l] = ch;
			l++;
		}
		else
		{
			if (l)
			{
				flag = false;
				word[l] = '\0';
			}
		}
		ipos++;
	} while (flag);
	var = atoi(word);
}

void get_double(char *str, int &ipos, double &var)
{
	int l;
	bool flag;
	char ch;
	char word[32];
	flag = true;
	l = 0;
	do
	{
		ch = str[ipos];
		if (!(ch == ' ' || ch == '\t' || ch == '\n'))
		{
			word[l] = ch;
			l++;
		}
		else
		{
			if (l)
			{
				flag = false;
				word[l] = '\0';
			}
		}
		ipos++;
	} while (flag);
	var = atof(word);
}

int main()
{
	int i,k,retp,m,m_end;
	double *B[2], *E[2], *A[2], *V[2], *H[2], *E2d[2];
	ifstream inf;
	ofstream ofp,ofps,ofpc;
	PointXYZ *xyzVectorE, *xyzVectorB;
	double x,y;
	FILE *fp;

	PointXYZ GenMin,GenMax;
	PointXYZ Av,Bv;

	int npntE,npntB;
	int *matrec;
	PointRZ *recE,*recB;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double nu,w;
	double DminSrc;

	int p1,p2,ipls,npls,nprof;
	vector<SqLoop> GenSq;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecToSourceB, RecToSourceE;

	RectOutData *vRecOutData;
	
	logfile.open("LogOutput2dEr");

	npls=0;

	inf.open("clcnplsa");
	if(!inf)
	{
		logfile<<"Error in open file "<<"clcnplsa"<<endl;
		cout<<"Error in open file "<<"clcnplsa"<<endl;
		return 1;
	}
	inf>>npls;
	inf.close();
	inf.clear();

	GenSq.resize(npls);

	inf.open("srsclca");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclca"<<endl;
		cout<<"Error in open file "<<"srsclca"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].A.x>>GenSq[i].A.y>>GenSq[i].A.z>>GenSq[i].B.x>>GenSq[i].B.y>>GenSq[i].B.z;}
	inf.close();
	inf.clear();

	RecvPlsIgB.resize(npls+1);
	RecvPlsIgE.resize(npls+1);

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

	for(i=0;i<npls;i++)
	{
		RecvPlsIgB[i+1]+=RecvPlsIgB[i];
		RecvPlsIgE[i+1]+=RecvPlsIgE[i];
	}

	npntB=RecvPlsIgB[npls];
	npntE=RecvPlsIgE[npls];

	RecToSourceB.resize(npntB);
	RecToSourceE.resize(npntE);

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

	task2d task("Ax");
	task2d taskout("Ax");
	SmoothOutput2D SobjB,SobjE;
	double rc0;

	retp=task.Read(1);
	if(retp)
	{
		logfile<<"Function Read() for Ax returned "<<retp<<'\n';
		return 1;
	}

	retp = taskout.Read(0);
	if(retp)
	{
		logfile << "Function Read() for AxOutput returned " << retp << '\n';
		return 1;
	}

	rc0 = 0.5*(taskout.rm[0] + taskout.rm[1]);

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
	for(i=0;i<npntB;i++)
	{
		Av=GenSq[RecToSourceB[i]].A;
		Bv=GenSq[RecToSourceB[i]].B;
		inf>>xyzVectorB[i].x>>xyzVectorB[i].y>>xyzVectorB[i].z;
		AddReciver(xyzVectorB[i],Av,recB[2*i],sinfiB[2*i],cosfiB[2*i],rc0);
		AddReciver(xyzVectorB[i],Bv,recB[2*i+1],sinfiB[2*i+1],cosfiB[2*i+1],rc0);
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
	xyzVectorE = new PointXYZ[npntE];
	recE = new PointRZ[npntE * 2];
	if(!recE)return RETCODE_NOMEM;
	sinfiE = new double[npntE * 2];
	if(!sinfiE)return RETCODE_NOMEM;
	cosfiE = new double[npntE * 2];
	if(!cosfiE)return RETCODE_NOMEM;
	matrec = new int[npntE * 2];
	if(!matrec)return RETCODE_NOMEM;	
	for (i = 0; i<npntE; i++)
	{
		Av=GenSq[RecToSourceE[i]].A;
		Bv=GenSq[RecToSourceE[i]].B;
		inf>>xyzVectorE[i].x>>xyzVectorE[i].y>>xyzVectorE[i].z;
		AddReciver(xyzVectorE[i],Av,recE[2*i],sinfiE[2*i],cosfiE[2*i],rc0);
		AddReciver(xyzVectorE[i],Bv,recE[2*i+1],sinfiE[2*i+1],cosfiE[2*i+1],rc0);
		matrec[2*i]=0;
		matrec[2*i+1]=0;
	}
	inf.close();
	inf.clear();

	A[0]=new double[task.kpnt*npls];
	if(!A[0])return RETCODE_NOMEM;
	fp=fopen("ar.re","rb");
	if(!fp) return RETCODE_NOFILE;
	fread(A[0],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	A[1]=new double[task.kpnt*npls];
	if(!A[1])return RETCODE_NOMEM;
	fp=fopen("ar.im","rb");
	if(!fp) return RETCODE_NOFILE;
	fread(A[1],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	V[0]=new double[task.kpnt*npls];
	if(!V[0])return RETCODE_NOMEM;
	fp=fopen("v.re","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(V[0],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	V[1]=new double[task.kpnt*npls];
	if(!V[1])return RETCODE_NOMEM;
	fp=fopen("v.im","rb");
	if(!fp)return RETCODE_NOFILE;
	fread(V[1],sizeof(double),task.kpnt*npls,fp);
	fclose(fp);

	B[0]=new double[npntB*2];
	if(!B[0])return RETCODE_NOMEM;
	B[1]=new double[npntB*2];
	if(!B[1])return RETCODE_NOMEM;

	E[0] = new double[npntE*2];
	if(!E[0])return RETCODE_NOMEM;
	E[1] = new double[npntE*2];
	if(!E[1])return RETCODE_NOMEM;

	H[0]=new double[task.kpnt*npls];
	if(!H[0])return RETCODE_NOMEM;
	H[1]=new double[task.kpnt*npls];
	if(!H[1])return RETCODE_NOMEM;

	for(i=0;i<task.kpnt*npls;i++)
	{
		H[0][i]=w*A[1][i];
		H[1][i]=-w*A[0][i];
	}

	retp=SobjB.Init(npntB*2,task.kpnt,task.krect,recB,task.pnt,task.rect,0,NULL,task.sigma,
		task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp<<endl;
		return retp;
	}
	
	SobjB.Output_dArdz(A[0],A[1],B[0],B[1],1,RecToSourceB);

	ofp.open("b2d_er");
	ofp << scientific << setprecision(14);
	for (i = 0; i<npntB; i++)
	{
		x = (-B[0][2 * i + 1] * sinfiB[2 * i + 1] + B[0][2 * i] * sinfiB[2 * i]);
		y = (B[0][2 * i + 1] * cosfiB[2 * i + 1] - B[0][2 * i] * cosfiB[2 * i]);
		ofp << x << ' ' << y << ' ' << 0.0 << ' ';
		x = (-B[1][2 * i + 1] * sinfiB[2 * i + 1] + B[1][2 * i] * sinfiB[2 * i]);
		y = (B[1][2 * i + 1] * cosfiB[2 * i + 1] - B[1][2 * i] * cosfiB[2 * i]);
		ofp << x << ' ' << y << ' ' << 0.0 << '\n';
	}
	ofp.close();
	ofp.clear();

	retp = SobjE.Init(npntE*2,task.kpnt,task.krect,recE,task.pnt,task.rect,0,NULL,task.sigma,
		task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm);
	if (retp)
	{
		logfile << "SobjE.Init() returned " << retp << endl;
		return retp;
	}

	SobjE.Output_Er(H[0],H[1],V[0],V[1],E[0],E[1],1,RecToSourceE);

	ofp.open("e2d_er");
	ofp << scientific << setprecision(14);
	for (i = 0; i<npntE; i++)
	{
		x = (E[0][2 * i + 1] * cosfiE[2 * i + 1] - E[0][2 * i] * cosfiE[2 * i]);
		y = (E[0][2 * i + 1] * sinfiE[2 * i + 1] - E[0][2 * i] * sinfiE[2 * i]);
		ofp << x << ' ' << y << ' ' << 0.0 << ' ';
		x = (E[1][2 * i + 1] * cosfiE[2 * i + 1] - E[1][2 * i] * cosfiE[2 * i]);
		y = (E[1][2 * i + 1] * sinfiE[2 * i + 1] - E[1][2 * i] * sinfiE[2 * i]);
		ofp << x << ' ' << y << ' ' << 0.0 << '\n';
	}
	ofp.close();
	ofp.clear();
	
	if (B[0]){ delete[] B[0]; B[0] = NULL; }
	if (B[1]){ delete[] B[1]; B[1] = NULL; }
	if (E[0]){ delete[] E[0]; E[0] = NULL; }
	if (E[1]){ delete[] E[1]; E[1] = NULL; }

	if (recE){ delete[] recE; recE = NULL; }
	if (recB){ delete[] recB; recB = NULL; }
	if (sinfiE){ delete[] sinfiE; sinfiE = NULL; }
	if (cosfiE){ delete[] cosfiE; cosfiE = NULL; }
	if (sinfiB){ delete[] sinfiB; sinfiB = NULL; }
	if (cosfiB){ delete[] cosfiB; cosfiB = NULL; }

	if (matrec){ delete[] matrec; matrec = NULL; }

	if (xyzVectorE){ delete[] xyzVectorE; xyzVectorE = NULL; }
	if (xyzVectorB){ delete[] xyzVectorB; xyzVectorB = NULL; }

	inf.open("xyzVectorE0", ios::binary);
	if(inf)
	{
		int ipos;
		char *sss;
		inf >> k;
		inf.seekg(0, inf.end);
		k = inf.tellg();
		inf.seekg(0, inf.beg);
		sss = new char[k];
		inf.read(sss, k);
		inf.close();

		ipos = 0;
		get_int(sss, ipos, npntE);

		xyzVectorE = new PointXYZ[npntE];
		recE = new PointRZ[npntE * 2];
		if (!recE)return RETCODE_NOMEM;
		sinfiE = new double[npntE * 2];
		if (!sinfiE)return RETCODE_NOMEM;
		cosfiE = new double[npntE * 2];
		if (!cosfiE)return RETCODE_NOMEM;
		matrec = new int[npntE * 2];
		if (!matrec)return RETCODE_NOMEM;

		for (i = 0; i < npntE; i++)
		{
			get_double(sss, ipos, xyzVectorE[i].x);
			get_double(sss, ipos, xyzVectorE[i].y);
			get_double(sss, ipos, xyzVectorE[i].z);
		}

		if (sss){ delete[] sss; sss = NULL; }
	}
	else
	{
		npntE=0;
	}
	inf.clear();

	inf.open("xyzVectorE0n");
	if(inf)
	{
		for (i = 0; i<npntE; i++)
		{
			inf>>matrec[2*i];
			matrec[2*i] = task.mtr3d2d[matrec[2*i]-1];
			matrec[2*i+1] = matrec[2*i];
		}
		inf.close();
	}
	inf.clear();

	E[0] = new double[npntE * 2];
	if (!E[0])return RETCODE_NOMEM;
	E[1] = new double[npntE * 2];
	if (!E[1])return RETCODE_NOMEM;

	E2d[0] = new double[taskout.kpnt * 2];
	if (!E2d[0])return RETCODE_NOMEM;
	E2d[1] = new double[taskout.kpnt * 2];
	if (!E2d[1])return RETCODE_NOMEM;

	vRecOutData = new RectOutData[npntE*2];
	if (!vRecOutData)return RETCODE_NOMEM;

	retp = SobjE.Init(taskout.qr, taskout.qz, task.kpnt, task.krect, taskout.pnt, task.pnt, task.rect, 0, NULL, task.sigma,
		task.nreg, task.qr, task.qz, task.reg, task.rm, task.zm);
	if(retp)
	{
		logfile << "SobjE.Init() returned " << retp << endl;
		return retp;
	}
	
	double rmin, rmax, zmin, zmax;

	if (npntE)
	{
		ofps.open("e_s_er.dat", ios::binary);
		ofpc.open("e_c_er.dat", ios::binary);

		for (ipls = 0; ipls < npls; ipls++)
		{
			SobjE.Output_Er(H[0] + task.kpnt*ipls, H[1] + task.kpnt*ipls, V[0] + task.kpnt*ipls, V[1] + task.kpnt*ipls, E2d[0], E2d[1], 1);

			for (i = 0; i<npntE; i++)
			{
				AddReciver(xyzVectorE[i], GenSq[ipls].A, recE[2 * i], sinfiE[2 * i], cosfiE[2 * i], rc0);
				AddReciver(xyzVectorE[i], GenSq[ipls].B, recE[2 * i + 1], sinfiE[2 * i + 1], cosfiE[2 * i + 1], rc0);
			}

			rmin = zmin = 1e+30;
			rmax = zmax = -1e+30;

			for (i = 0; i < npntE * 2; i++)
			{
				if (recE[i].r < rmin){ rmin = recE[i].r; }
				if (recE[i].r > rmax){ rmax = recE[i].r; }
				if (recE[i].z < zmin){ zmin = recE[i].z; }
				if (recE[i].z > zmax){ zmax = recE[i].z; }
			}

			retp = FindElementsForReceivers(npntE*2, taskout.qr, taskout.qz, taskout.rm, taskout.zm, taskout.reg, recE, vRecOutData, taskout.pnt, taskout.rect, rmin, zmin, rmax, zmax);
			if (retp)
			{
				logfile << "FindElementsForReceivers returned " << retp << endl;
				return retp;
			}

			Output_Field(npntE*2, E2d[0], E2d[1], E[0], E[1], vRecOutData);


			for (i = 0; i<npntE; i++)
			{
				x = (E[0][2 * i + 1] * cosfiE[2 * i + 1] - E[0][2 * i] * cosfiE[2 * i]);
				y = (E[0][2 * i + 1] * sinfiE[2 * i + 1] - E[0][2 * i] * sinfiE[2 * i]);
				if (recE[i].r>DminSrc || (fabs(recE[i].z - GenSq[ipls].A.z)>DminSrc && fabs(recE[i].z - GenSq[ipls].B.z)>DminSrc))
				{
					ofps<x<y<0.0;
				}
				else
				{
					ofps<0.0<0.0<0.0;
				}
				x = (E[1][2 * i + 1] * cosfiE[2 * i + 1] - E[1][2 * i] * cosfiE[2 * i]);
				y = (E[1][2 * i + 1] * sinfiE[2 * i + 1] - E[1][2 * i] * sinfiE[2 * i]);
				if (recE[i].r>DminSrc || (fabs(recE[i].z - GenSq[ipls].A.z)>DminSrc && fabs(recE[i].z - GenSq[ipls].B.z)>DminSrc))
				{
					ofpc<x<y<0.0;
				}
				else
				{
					ofpc<0.0<0.0<0.0;
				}
			}
		}

		ofps.close();
		ofps.clear();
		ofpc.close();
		ofpc.clear();
	}

	if (E2d[0]){ delete[] E2d[0]; E2d[0] = NULL; }
	if (E2d[1]){ delete[] E2d[1]; E2d[1] = NULL; }

	if(A[0]){delete [] A[0]; A[0]=NULL;}
	if(A[1]){delete [] A[1]; A[1]=NULL;}
	if(E[0]){delete [] E[0]; E[0]=NULL;}
	if(E[1]){delete [] E[1]; E[1]=NULL;}
	if(V[0]){delete [] V[0]; V[0]=NULL;}
	if(V[1]){delete [] V[1]; V[1]=NULL;}
	if(H[0]){delete [] H[0]; H[0]=NULL;}
	if(H[1]){delete [] H[1]; H[1]=NULL;}

	if(recE){delete [] recE; recE=NULL;}
	if(sinfiE){delete [] sinfiE; sinfiE=NULL;}
	if(cosfiE){delete [] cosfiE; cosfiE=NULL;}

	if (matrec){ delete[] matrec; matrec = NULL; }
	if (vRecOutData){ delete[] vRecOutData; vRecOutData = NULL; }
	if (xyzVectorE){ delete[] xyzVectorE; xyzVectorE = NULL; }

	logfile.close();
	logfile.clear();

	return 0;
}
