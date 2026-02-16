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

struct Line
{
	PointXYZ A,B;
};

ofstream logfile;

void GetSource(PointXYZ &GenMin,PointXYZ &GenMax,PointXYZ &pA,PointXYZ &pB,PointXYZ &vAB,PointXYZ &Ic,double &len,double &cosphi,double &sinphi)
{
	double phi,cedr;
	PointXYZ pT;

	pA.x=GenMin.x;
	pA.y=GenMin.y;
	pA.z=GenMin.z;

	pB.x=GenMax.x;
	pB.y=GenMax.y;
	pB.z=GenMin.z;

	vAB.x=pB.x-pA.x;
	vAB.y=pB.y-pA.y;
	vAB.z=pB.z-pA.z;

	cedr=sqrt(vAB.x*vAB.x+vAB.y*vAB.y+vAB.z*vAB.z);
	phi=acos(vAB.x/cedr);

	if(vAB.y<0) phi*=-1.;

	cosphi=cos(-phi);
	sinphi=sin(-phi);
	
	pB.x=pB.x-pA.x;
	pB.y=pB.y-pA.y;

	pT.x=pB.x*cosphi-pB.y*sinphi;
	pT.y=pB.x*sinphi+pB.y*cosphi;

	pB.x=pT.x+pA.x;
	pB.y=pT.y+pA.y;

	vAB.x=pB.x-pA.x;
	vAB.y=pB.y-pA.y;
	vAB.z=pB.z-pA.z;

	Ic=vAB;

	Ic.x/=cedr;
	Ic.y/=cedr;
	Ic.z/=cedr;

	len=cedr/NUMBEROFABDIPOLES;

	vAB.x/=NUMBEROFABDIPOLES;
	vAB.y/=NUMBEROFABDIPOLES;
	vAB.z/=NUMBEROFABDIPOLES;
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
	int i,idp,k,retp,m,m_end;
	double *Bx[2],*By[2],*Bz[2],*E[2],*E2d[2];
	double U[3][2],ux,uy,uz;
	ifstream inf;
	ofstream ofp,ofps,ofpc;
	PointXYZ *xyzVectorE,*xyzVectorB;
	double x,y,r,len;
	double cosphi,sinphi;
	double nu,w;

	int npntE,npntB,npntE0;
	int *matrec;
	PointRZ *recE,*recB;
	double *sinfiE,*cosfiE;
	double *sinfiB,*cosfiB;
	double DminSrc;
	
	int p1,p2,ipls,npls,nprof;
	vector<Line> GenLin;
	vector<int> RecvPlsIgB, RecvPlsIgE, RecToSourceB, RecToSourceE;
	vector<double> cosphiB,sinphiB,cosphiE,sinphiE,lenE,lenB;
	
	RectOutData *vRecOutData;

	NUMBEROFABDIPOLES=20;

	logfile.open("LogOutput2dAx");

	inf.open("numberofabdipoles");
	if(!inf)
	{
		logfile<<"Error in open file "<<"numberofabdipoles"<<endl;
		cout<<"Error in open file "<<"numberofabdipoles"<<endl;
		return 1;
	}
	inf>>NUMBEROFABDIPOLES;
	inf.close();
	inf.clear();

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

	GenLin.resize(npls);

	inf.open("srsclca");
	if(!inf)
	{
		logfile<<"Error in open file "<<"srsclca"<<endl;
		cout<<"Error in open file "<<"srsclca"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenLin[i].A.x>>GenLin[i].A.y>>GenLin[i].A.z>>GenLin[i].B.x>>GenLin[i].B.y>>GenLin[i].B.z;}
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

	PointXYZ pA,pB,pT,vAB,Ic;
	vector<PointXYZ> IcE,IcB;

	IcB.resize(npntB);
	IcE.resize(npntE);
	cosphiB.resize(npntB);
	sinphiB.resize(npntB);
	cosphiE.resize(npntE);
	sinphiE.resize(npntE);
	lenB.resize(npntB);
	lenE.resize(npntE);

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

	retp = task.ReadSolution(npls);
	if (retp)
	{
		logfile << "Function Read() for Ax returned " << retp << '\n';
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
	recB=new PointRZ[npntB*NUMBEROFABDIPOLES];
	if(!recB)return RETCODE_NOMEM;
	sinfiB=new double[npntB*NUMBEROFABDIPOLES];
	if(!sinfiB)return RETCODE_NOMEM;
	cosfiB=new double[npntB*NUMBEROFABDIPOLES];
	if(!cosfiB)return RETCODE_NOMEM;
	for(i=0;i<npntB;i++)
	{
		GetSource(GenLin[RecToSourceB[i]].A,GenLin[RecToSourceB[i]].B,pA,pB,vAB,Ic,len,cosphi,sinphi);
		IcB[i]=Ic;
		cosphiB[i]=cosphi;
		sinphiB[i]=sinphi;
		lenB[i]=len;
		inf>>xyzVectorB[i].x>>xyzVectorB[i].y>>xyzVectorB[i].z;
		for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
		{
			k=i*NUMBEROFABDIPOLES+idp;
			x=xyzVectorB[i].x;
			y=xyzVectorB[i].y;

			x=x-pA.x;
			y=y-pA.y;
			pT.x=x*cosphi-y*sinphi;
			pT.y=x*sinphi+y*cosphi;	
			x=pT.x+pA.x;
			y=pT.y+pA.y;

			pT.x=pA.x+vAB.x*(idp+0.5);
			pT.y=pA.y+vAB.y*(idp+0.5);
			x-=pT.x;
			y-=pT.y;
			r=sqrt(x*x+y*y);
			if(r<rc0)r=rc0;
			recB[k].r=r;
			recB[k].z=xyzVectorB[i].z;
			sinfiB[k]=y/r;
			cosfiB[k]=x/r;
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
	recE=new PointRZ[npntE*NUMBEROFABDIPOLES];
	if(!recE)return RETCODE_NOMEM;
	sinfiE=new double[npntE*NUMBEROFABDIPOLES];
	if(!sinfiE)return RETCODE_NOMEM;
	cosfiE=new double[npntE*NUMBEROFABDIPOLES];
	if(!cosfiE)return RETCODE_NOMEM;
	matrec=new int[npntE*NUMBEROFABDIPOLES];
	if(!matrec)return RETCODE_NOMEM;
	for (i = 0; i<npntE; i++)
	{
		GetSource(GenLin[RecToSourceE[i]].A,GenLin[RecToSourceE[i]].B,pA,pB,vAB,Ic,len,cosphi,sinphi);
		IcE[i]=Ic;
		cosphiE[i]=cosphi;
		sinphiE[i]=sinphi;
		lenE[i]=len;
		inf>>xyzVectorE[i].x>>xyzVectorE[i].y>>xyzVectorE[i].z;
		for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
		{
			k=i*NUMBEROFABDIPOLES+idp;
			x=xyzVectorE[i].x;
			y=xyzVectorE[i].y;

			x=x-pA.x;
			y=y-pA.y;
			pT.x=x*cosphi-y*sinphi;
			pT.y=x*sinphi+y*cosphi;	
			x=pT.x+pA.x;
			y=pT.y+pA.y;

			pT.x=pA.x+vAB.x*(idp+0.5);
			pT.y=pA.y+vAB.y*(idp+0.5);
			x-=pT.x;
			y-=pT.y;
			r=sqrt(x*x+y*y);
			if(r<rc0)r=rc0;
			recE[k].r=r;
			recE[k].z=xyzVectorE[i].z;
			sinfiE[k]=y/r;
			cosfiE[k]=x/r;
			matrec[k]=0;
		}
	}
	inf.close();
	inf.clear();

	Bx[0]=new double[npntB*NUMBEROFABDIPOLES];
	if(!Bx[0])return RETCODE_NOMEM;
	Bx[1]=new double[npntB*NUMBEROFABDIPOLES];
	if(!Bx[1])return RETCODE_NOMEM;

	By[0]=new double[npntB*NUMBEROFABDIPOLES];
	if(!By[0])return RETCODE_NOMEM;
	By[1]=new double[npntB*NUMBEROFABDIPOLES];
	if(!By[1])return RETCODE_NOMEM;

	Bz[0]=new double[npntB*NUMBEROFABDIPOLES];
	if(!Bz[0])return RETCODE_NOMEM;
	Bz[1]=new double[npntB*NUMBEROFABDIPOLES];
	if(!Bz[1])return RETCODE_NOMEM;

	E[0] = new double[npntE*NUMBEROFABDIPOLES];
	if(!E[0])return RETCODE_NOMEM;
	E[1] = new double[npntE*NUMBEROFABDIPOLES];
	if(!E[1])return RETCODE_NOMEM;

	retp=SobjB.Init(npntB*NUMBEROFABDIPOLES,task.kpnt,task.krect,recB,task.pnt,task.rect,0,NULL,task.sigma,
		task.nreg,task.qr,task.qz,task.reg,task.rm,task.zm);
	if(retp)
	{
		logfile<<"SobjB.Init() returned "<<retp;
		return retp;
	}
	SobjB.OutputB(task.v2s,task.v2c,Bx[0],Bx[1],By[0],By[1],Bz[0],Bz[1],IcB,sinfiB,cosfiB,1,RecToSourceB);

	ofp.open("b2d_ax");
	ofp<<scientific<<setprecision(14);
	for(i=0;i<npntB;i++)
	{
		Ic=IcB[i];
		cosphi=cosphiB[i];
		sinphi=-sinphiB[i];
		len=lenB[i];

		U[0][0]=U[1][0]=U[2][0]=0.0;
		U[0][1]=U[1][1]=U[2][1]=0.0;
		for(idp=0;idp<NUMBEROFABDIPOLES;idp++)
		{
			k=i*NUMBEROFABDIPOLES+idp;
			U[0][0]+=Bx[0][k];
			U[0][1]+=Bx[1][k];
			U[1][0]+=By[0][k];
			U[1][1]+=By[1][k];
			U[2][0]+=Bz[0][k];
			U[2][1]+=Bz[1][k];
		}
		U[0][0]*=len;
		U[0][1]*=len;
		U[1][0]*=len;
		U[1][1]*=len;
		U[2][0]*=len;
		U[2][1]*=len;

		ux=U[0][0];
		uy=U[1][0];
		uz=U[2][0];

		x=ux*cosphi-uy*sinphi;
		y=ux*sinphi+uy*cosphi;

		ux=x;
		uy=y;

		ofp<<ux<<' '<<uy<<' '<<uz<<' ';

		ux=U[0][1];
		uy=U[1][1];
		uz=U[2][1];

		x=ux*cosphi-uy*sinphi;
		y=ux*sinphi+uy*cosphi;

		ux=x;
		uy=y;

		ofp<<ux<<' '<<uy<<' '<<uz<<'\n';
	}
	ofp.close();
	ofp.clear();
		
	retp = SobjE.Init(npntE*NUMBEROFABDIPOLES, task.kpnt, task.krect, recE, task.pnt, task.rect, 0, NULL, task.sigma,
		task.nreg, task.qr, task.qz, task.reg, task.rm, task.zm);
	if (retp)
	{
		logfile << "SobjE.Init() returned " << retp;
		return retp;
	}
	SobjE.OutputE(task.v2s, task.v2c, E[0], E[1], 1, RecToSourceE);

	ofp.open("e2d_ax");
	ofp << scientific << setprecision(14);

	for (i = 0; i<npntE; i++)
	{
		Ic = IcE[i];
		cosphi = cosphiE[i];
		sinphi = -sinphiE[i];
		len = lenE[i];

		U[0][0] = U[0][1] = 0.0;
		for (idp = 0; idp<NUMBEROFABDIPOLES; idp++)
		{
			k = i*NUMBEROFABDIPOLES + idp;
			U[0][0] += E[0][k];
			U[0][1] += E[1][k];
		}
		U[0][0] *= len;
		U[0][1] *= len;

		ux = U[0][1] * Ic.x;
		uy = U[0][1] * Ic.y;
		uz = U[0][1] * Ic.z;

		x = ux*cosphi - uy*sinphi;
		y = ux*sinphi + uy*cosphi;

		ux = x;
		uy = y;

		ofp << w*ux << ' ' << w*uy << ' ' << w*uz << ' ';

		ux = U[0][0] * Ic.x;
		uy = U[0][0] * Ic.y;
		uz = U[0][0] * Ic.z;

		x = ux*cosphi - uy*sinphi;
		y = ux*sinphi + uy*cosphi;

		ux = x;
		uy = y;

		ofp << -w*ux << ' ' << -w*uy << ' ' << -w*uz << '\n';
	}
	ofp.close();
	ofp.clear();

	if (Bx[0]){ delete[] Bx[0]; Bx[0] = NULL; }
	if (Bx[1]){ delete[] Bx[1]; Bx[1] = NULL; }
	if (By[0]){ delete[] By[0]; By[0] = NULL; }
	if (By[1]){ delete[] By[1]; By[1] = NULL; }
	if (Bz[0]){ delete[] Bz[0]; Bz[0] = NULL; }
	if (Bz[1]){ delete[] Bz[1]; Bz[1] = NULL; }

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
	if (inf)
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
		recE = new PointRZ[npntE*NUMBEROFABDIPOLES];
		if (!recE)return RETCODE_NOMEM;
		sinfiE = new double[npntE*NUMBEROFABDIPOLES];
		if (!sinfiE)return RETCODE_NOMEM;
		cosfiE = new double[npntE*NUMBEROFABDIPOLES];
		if (!cosfiE)return RETCODE_NOMEM;
		matrec = new int[npntE*NUMBEROFABDIPOLES];
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
		npntE = 0;
	}
	inf.clear();

	inf.open("xyzVectorE0n");
	if (inf)
	{
		for (i = 0; i<npntE; i++)
		{
			inf >> matrec[i*NUMBEROFABDIPOLES];
			for (idp = 0; idp<NUMBEROFABDIPOLES; idp++)
			{
				k = i*NUMBEROFABDIPOLES + idp;
				matrec[k] = task.mtr3d2d[matrec[i*NUMBEROFABDIPOLES] - 1];
			}
		}
		inf.close();
	}
	inf.clear();

	IcE.resize(npntE);
	cosphiE.resize(npntE);
	sinphiE.resize(npntE);
	lenE.resize(npntE);

	E[0] = new double[npntE * NUMBEROFABDIPOLES];
	if (!E[0])return RETCODE_NOMEM;
	E[1] = new double[npntE * NUMBEROFABDIPOLES];
	if (!E[1])return RETCODE_NOMEM;

	E2d[0] = new double[taskout.kpnt * 2];
	if (!E2d[0])return RETCODE_NOMEM;
	E2d[1] = new double[taskout.kpnt * 2];
	if (!E2d[1])return RETCODE_NOMEM;

	vRecOutData = new RectOutData[npntE * NUMBEROFABDIPOLES];
	if (!vRecOutData)return RETCODE_NOMEM;
	
	retp = SobjE.Init(taskout.qr, taskout.qz, task.kpnt, task.krect, taskout.pnt, task.pnt, task.rect, 0, NULL, task.sigma,
		task.nreg, task.qr, task.qz, task.reg, task.rm, task.zm);
	if (retp)
	{
		logfile << "SobjE.Init() returned " << retp << endl;
		return retp;
	}

	double rmin, rmax, zmin, zmax;

	if (npntE)
	{
		ofps.open("e_s_ax.dat", ios::binary);
		ofpc.open("e_c_ax.dat", ios::binary);

		for (ipls = 0; ipls < npls; ipls++)
		{
			SobjE.OutputE(task.v2s + task.kpnt*ipls, task.v2c + task.kpnt*ipls, E2d[0], E2d[1], 1);

			for (i = 0; i<npntE; i++)
			{

				GetSource(GenLin[ipls].A, GenLin[ipls].B, pA, pB, vAB, Ic, len, cosphi, sinphi);
				IcE[i] = Ic;
				cosphiE[i] = cosphi;
				sinphiE[i] = sinphi;
				lenE[i] = len;
				for (idp = 0; idp<NUMBEROFABDIPOLES; idp++)
				{
					k = i*NUMBEROFABDIPOLES + idp;
					x = xyzVectorE[i].x;
					y = xyzVectorE[i].y;

					x = x - pA.x;
					y = y - pA.y;
					pT.x = x*cosphi - y*sinphi;
					pT.y = x*sinphi + y*cosphi;
					x = pT.x + pA.x;
					y = pT.y + pA.y;

					pT.x = pA.x + vAB.x*(idp + 0.5);
					pT.y = pA.y + vAB.y*(idp + 0.5);
					x -= pT.x;
					y -= pT.y;
					r = sqrt(x*x + y*y);
					if (r<rc0)r = rc0;
					recE[k].r = r;
					recE[k].z = xyzVectorE[i].z;
					sinfiE[k] = y / r;
					cosfiE[k] = x / r;
					matrec[k] = 0;
				}
			}

			rmin = zmin = 1e+30;
			rmax = zmax = -1e+30;

			for (i = 0; i < npntE * NUMBEROFABDIPOLES; i++)
			{
				if (recE[i].r < rmin){ rmin = recE[i].r; }
				if (recE[i].r > rmax){ rmax = recE[i].r; }
				if (recE[i].z < zmin){ zmin = recE[i].z; }
				if (recE[i].z > zmax){ zmax = recE[i].z; }
			}

			retp = FindElementsForReceivers(npntE*NUMBEROFABDIPOLES, taskout.qr, taskout.qz, taskout.rm, taskout.zm, taskout.reg, recE, vRecOutData, taskout.pnt, taskout.rect, rmin, zmin, rmax, zmax);
			if (retp)
			{
				logfile << "FindElementsForReceivers returned " << retp << endl;
				return retp;
			}
			Output_Field(npntE*NUMBEROFABDIPOLES, E2d[0], E2d[1], E[0], E[1], vRecOutData);

			for (i = 0; i<npntE; i++)
			{
				Ic = IcE[i];
				cosphi = cosphiE[i];
				sinphi = -sinphiE[i];
				len = lenE[i];

				U[0][0] = U[0][1] = 0.0;
				for (idp = 0; idp<NUMBEROFABDIPOLES; idp++)
				{
					k = i*NUMBEROFABDIPOLES + idp;
					U[0][0] += E[0][k];
					U[0][1] += E[1][k];
				}
				U[0][0] *= len;
				U[0][1] *= len;

				ux = U[0][1] * Ic.x;
				uy = U[0][1] * Ic.y;
				uz = U[0][1] * Ic.z;

				x = ux*cosphi - uy*sinphi;
				y = ux*sinphi + uy*cosphi;

				ux = x;
				uy = y;

				if (recE[i].r>DminSrc || (fabs(recE[i].z - GenLin[ipls].A.z)>DminSrc && fabs(recE[i].z - GenLin[ipls].B.z)>DminSrc))
				{
					ofps<w*ux<w*uy<w*uz;
				}
				else
				{
					ofps<0.0<0.0<0.0;
				}

				ux = U[0][0] * Ic.x;
				uy = U[0][0] * Ic.y;
				uz = U[0][0] * Ic.z;

				x = ux*cosphi - uy*sinphi;
				y = ux*sinphi + uy*cosphi;

				ux = x;
				uy = y;

				if (recE[i].r>DminSrc || (fabs(recE[i].z - GenLin[ipls].A.z)>DminSrc && fabs(recE[i].z - GenLin[ipls].B.z)>DminSrc))
				{
					ofpc<-w*ux<-w*uy<-w*uz;
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

	if(E[0]){delete [] E[0]; E[0]=NULL;}
	if(E[1]){delete [] E[1]; E[1]=NULL;}
	
	if(recE){delete [] recE; recE=NULL;}
	if(sinfiE){delete [] sinfiE; sinfiE=NULL;}
	if(cosfiE){delete [] cosfiE; cosfiE=NULL;}

	if(matrec){delete [] matrec; matrec = NULL; }

	if(vRecOutData){delete [] vRecOutData; vRecOutData = NULL; }

	if(xyzVectorE){delete [] xyzVectorE; xyzVectorE = NULL; }

	logfile.close();
	logfile.clear();

	return 0;
}
