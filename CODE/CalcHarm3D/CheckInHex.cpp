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
 *  This file contains code of the functions to test for hit inside the hexahedron
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024   
*/

#include "stdafx.h"
#include "CheckInHex.h"

extern ofstream logfile;

int slae(double m[][3],double f[])
{
	int i,j;
	double det,im[3][3],x[3],z;

	det=GetDeterminant33(m);

	if(fabs(det)>1e-6)
	{
		InverseMatrix33WithDet(m,im,det);
	}
	else
	{
		logfile << "Regularization Jacobian matrix" << '\n';

		det=0.0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				z=fabs(m[i][j]);
				det=(det>z)? det : z;
			}
		}

		det=(det>1e-6)? det : 1e-6;

		for(i=0;i<3;i++)
		{
			m[i][i]+=1e-6*det;
		}

		InverseMatrix33(m,im);
	}

	for(i=0;i<3;i++)
	{
		x[i]=f[i];
		f[i]=0.0;
	}

	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			f[i]+=im[i][j]*x[j];
		}
	}

	return 0;
}

int gauss(int n,double *M,double *F)
{
	int ki,kj,kl,kin,kln;
	double aii,ali,fi;
	for(ki=0;ki<n;ki++)
	{
		kin=ki*n;
		aii=M[kin+ki];
		if(aii>-1e-20 && aii<1e-20)
		{
			logfile<<"Regularization Jacobian matrix"<<'\n';
			aii=1e-20;
		}
		aii=1.0/aii;
		M[kin+ki]=1.0;
		F[ki]*=aii;
		fi=F[ki];
		for(kj=n-1;kj>ki;kj--){M[kin+kj]*=aii;}
		for(kl=ki+1;kl<n;kl++)
		{
			kln=kl*n;
			ali=M[kln+ki];
			F[kl]-=ali*fi;
			for(kj=n-1;kj>ki;kj--){M[kln+kj]-=ali*M[kin+kj];}
		}
	}
	for(ki=n-1;ki>=0;ki--)
	{
		fi=F[ki];
		for(kl=ki-1;kl>=0;kl--){F[kl]-=M[kl*n+ki]*fi;}
	}
	return 0;
}

int inverse3(double m[3][3], double r[3][3], double epsdet)
{
	int i, j, k, li, lj, ii, jj;
	double a[2][2], det, sig, sigi, sigj;
	det = 0.0;
	sig = 1.0;

	for (k = 0; k < 3; k++)
	{
		for (lj = j = 0; j < 3; j++)
		{
			if (j != k)
			{
				for (i = 1; i < 3; i++)
				{
					a[i - 1][lj] = m[i][j];
				}
				lj++;
			}
		}
		det += sig * (a[0][0] * a[1][1] - a[0][1] * a[1][0]) * m[0][k];
		sig = -sig;
	}

	if (det > -epsdet && det < epsdet)
	{
		logfile << "det= " << det << '\n';
		return 1;
	}

	sigi = 1.0;
	for (ii = 0; ii < 3; ii++)
	{
		sigj = sigi;
		for (jj = 0; jj < 3; jj++)
		{
			for (li = i = 0; i < 3; i++)
			{
				if (i != ii)
				{
					for (lj = j = 0; j < 3; j++)
					{
						if (j != jj)
						{
							a[li][lj] = m[i][j];
							lj++;
						}
					}
					li++;
				}
			}

			r[jj][ii] = sigj * (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
			sigj = -sigj;
		}
		sigi = -sigi;
	}

	return 0;
}

int signum(double x){return ((x>0.0)? 1 : (x<0.0)? -1 : 0);}

double l0(double x){return 1.0-x;}
double l1(double x){return x;}

double phi(int i,double x,double y,double z)
{
	double val=0.0;
	switch(i)
	{
		case 0 : {val=l0(x)*l0(y)*l0(z);break;}
		case 1 : {val=l1(x)*l0(y)*l0(z);break;}
		case 2 : {val=l0(x)*l1(y)*l0(z);break;}
		case 3 : {val=l1(x)*l1(y)*l0(z);break;}
		case 4 : {val=l0(x)*l0(y)*l1(z);break;}
		case 5 : {val=l1(x)*l0(y)*l1(z);break;}
		case 6 : {val=l0(x)*l1(y)*l1(z);break;}
		case 7 : {val=l1(x)*l1(y)*l1(z);break;}
	}
	return val;
}

double dphidx(int i,double x,double y,double z)
{
	double val=0.0;
	switch(i)
	{
		case 0 : {val=-l0(y)*l0(z);break;}
		case 1 : {val=l0(y)*l0(z);break;}
		case 2 : {val=-l1(y)*l0(z);break;}
		case 3 : {val=l1(y)*l0(z);break;}
		case 4 : {val=-l0(y)*l1(z);break;}
		case 5 : {val=l0(y)*l1(z);break;}
		case 6 : {val=-l1(y)*l1(z);break;}
		case 7 : {val=l1(y)*l1(z);break;}
	}
	return val;
}

double dphidy(int i,double x,double y,double z)
{
	double val=0.0;
	switch(i)
	{
		case 0 : {val=-l0(x)*l0(z);break;}
		case 1 : {val=-l1(x)*l0(z);break;}
		case 2 : {val=l0(x)*l0(z);break;}
		case 3 : {val=l1(x)*l0(z);break;}
		case 4 : {val=-l0(x)*l1(z);break;}
		case 5 : {val=-l1(x)*l1(z);break;}
		case 6 : {val=l0(x)*l1(z);break;}
		case 7 : {val=l1(x)*l1(z);break;}
	}
	return val;
}

double dphidz(int i,double x,double y,double z)
{
	double val=0.0;
	switch(i)
	{
		case 0 : {val=-l0(x)*l0(y);break;}
		case 1 : {val=-l1(x)*l0(y);break;}
		case 2 : {val=-l0(x)*l1(y);break;}
		case 3 : {val=-l1(x)*l1(y);break;}
		case 4 : {val=l0(x)*l0(y);break;}
		case 5 : {val=l1(x)*l0(y);break;}
		case 6 : {val=l0(x)*l1(y);break;}
		case 7 : {val=l1(x)*l1(y);break;}
	}
	return val;
}

void HexMap(Point3D &lp,Point3D *Hex,Point3D &gp)
{
	int i;
	double val;
	gp.x=0.0;
	gp.y=0.0;
	gp.z=0.0;
	for(i=0;i<8;i++)
	{
		val=phi(i,lp.x,lp.y,lp.z);
		gp.x+=val*Hex[i].x;
		gp.y+=val*Hex[i].y;
		gp.z+=val*Hex[i].z;
	}
}

void CalcJ(Point3D &lp,Point3D *Hex,Point3D &gp,double J[3][3])
{
	double mdphidc[3][8],crd[3][8],cur[3],glb[3],a,b;
	int i,j,k;
	Point3D cp;
	for(k=0;k<8;k++)
	{
		crd[0][k]=Hex[k].x;
		crd[1][k]=Hex[k].y;
		crd[2][k]=Hex[k].z;
		mdphidc[0][k]=dphidx(k,lp.x,lp.y,lp.z);
		mdphidc[1][k]=dphidy(k,lp.x,lp.y,lp.z);
		mdphidc[2][k]=dphidz(k,lp.x,lp.y,lp.z);
	}
	HexMap(lp,Hex,cp);
	cur[0]=cp.x;	cur[1]=cp.y;	cur[2]=cp.z;
	glb[0]=gp.x;	glb[1]=gp.y;	glb[2]=gp.z;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			b=0.0;
			for(k=0;k<8;k++){b+=mdphidc[j][k]*crd[i][k];}
			a=signum(cur[i]-glb[i]);
			J[i][j]=a*b;
		}
	}
}

void CalcF(Point3D &lp,Point3D *Hex,Point3D &gp,double *F)
{
	Point3D cp;
	HexMap(lp,Hex,cp);
	F[0]=fabs(cp.x-gp.x);
	F[1]=fabs(cp.y-gp.y);
	F[2]=fabs(cp.z-gp.z);
}

double CalcVectorNorm(int n,double *F)
{
	int i;
	double res;
	res=0.0;
	for(i=0;i<n;i++){res+=F[i]*F[i];}
	return res;
}

int CheckInHex(Point3D *Hex,Point3D &gp,Point3D &lp)
{
	Point3D cp;
	double J[3][3],F[3],x[3],xo[3],Jdi[3],S[3], Jinv[3][3];

	int i,j,n,iter,maxiter;
	double res,lam,eps,resc;

	n=3;

	eps=1e-6;

	iter=0;
	maxiter=20;

	CalcF(lp,Hex,gp,F);
	res=CalcVectorNorm(n,F);
	x[0]=lp.x;
	x[1]=lp.y;
	x[2]=lp.z;

	while(res>eps)
	{
		for(i=0;i<n;i++)xo[i]=x[i];
		CalcJ(lp,Hex,gp,J);
		for(i=0;i<n;i++){Jdi[i]=J[i][i];}
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (fabs(J[i][j]) > 1e-12)
				{
					break;
				}
			}
			if (j == n)
			{
				F[i] = 0.0;
				for (j = 0; j < n; j++)
				{
					J[i][j] = (j != i) ? 0.0 : 1.0;
				}
			}
		}
		if(!gauss(n, &(J[0][0]), F) && iter<maxiter/2)
		{
			for(i=0;i<n;i++)x[i]=xo[i]-F[i];
			lp.x=x[0];
			lp.y=x[1];
			lp.z=x[2];
			CalcF(lp,Hex,gp,F);
			resc=CalcVectorNorm(n,F);
		}
		else
		{
			lam=1e-6;
			do{
				for(i=0;i<n;i++){x[i]=xo[i]-lam*Jdi[i];}
				lp.x=x[0];
				lp.y=x[1];
				lp.z=x[2];
				CalcF(lp,Hex,gp,F);
				resc=CalcVectorNorm(n,F);
				lam*=0.5;
			}while(!(resc<res) && lam>1e-14);
		}

		res=resc;

		iter++;
		if(iter>maxiter)
		{
			break;
		}

		for(i=0;i<n;i++)S[i]=x[i]-xo[i];
		if(CalcVectorNorm(n,S)<1e-6){break;}
	}

	HexMap(lp,Hex,cp);

	if (res > 1e-4)
	{
		logfile << '\n'
			<< res << " " << iter << " " << maxiter << '\n'
			<< gp.x << " " << gp.y << " " << gp.z << '\n'
			<< cp.x << " " << cp.y << " " << cp.z << '\n'
			<< lp.x << " " << lp.y << " " << lp.z << '\n'
			<< '\n';
	}

	return (iter>maxiter);
}
