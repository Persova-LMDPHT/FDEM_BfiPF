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
 *  This file contains code for smooth FEM output E and B fields in receivers.
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

extern ofstream logfile;

SmoothOutput2D::SmoothOutput2D()
{
	nrec=kpnt=krect=0;
	rec=NULL;
	pnt=NULL;
	rect=NULL;
	OutRecList=NULL;
	RecToElem.clear();
}

SmoothOutput2D::~SmoothOutput2D()
{
	RecToElem.clear();
	if(OutRecList){delete [] OutRecList; OutRecList=NULL;}
}

// Initialization of smoothing output (find elements for receivers)
int SmoothOutput2D::Init(int rec_qr, int rec_qz,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
						 int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm)
{
	int i, j, k, l, m, ii, jj, ef,iii,jjj;
	PointRZ celm;

	bool fstop;

	nrec = rec_qr*rec_qz;
	kpnt=p_kpnt;
	krect=p_krect;
	rec=p_rec;
	pnt=p_pnt;
	rect=p_rect;
	sigma=p_sigma;
	nreg=p_nreg;
	qr=p_qr;
	qz=p_qz;
	reg=p_reg;
	rm=p_rm;
	zm=p_zm;

	if(nrec==0)return 0;

	if ((int)RecToElem.size())RecToElem.clear();
	if (OutRecList){delete[] OutRecList; OutRecList = NULL;}

	RecToElem.resize(nrec);
	OutRecList=new Output2D[nrec];

	fstop=false;
	for (jjj=0;jjj<rec_qz;jjj++)
	{
		for(iii=0;iii<rec_qr;iii++)
		{
			k = jjj*rec_qr + iii;

			i = iii;
			i -= (i>0);
			j = jjj;
			j -= (j>0);

			m=j*(qr-1)+i;
			l=reg[m]-1;
			RecToElem[k]=l;

			celm.r=0.0;
			celm.z=0.0;
			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].oror.lem[ii]=-1;
				OutRecList[k].oz.lem[ii]=-1;
				jj=rect[l].nodes[ii]-1;
				celm.r+=pnt[jj].r;
				celm.z+=pnt[jj].z;
			}
			celm.r*=0.25;
			celm.z*=0.25;

			ef=!(rec[k].r>celm.r);
			GetSmoothForR(i,j,k,ef,l,0,rec[k],OutRecList[k].oror.lem,OutRecList[k].oror.coef,OutRecList[k].oror.v,OutRecList[k].oror.ves,OutRecList[k].oror.fl,0);

			ef=!(rec[k].z>celm.z);
			GetSmoothForZ(i,j,k,ef,l,0,rec[k],OutRecList[k].oz.lem,OutRecList[k].oz.coef,OutRecList[k].oz.v,OutRecList[k].oz.ves,OutRecList[k].oz.fl,0);
		}
	}

	return fstop;
}
int SmoothOutput2D::Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
						 int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm)
{
	int i,j,k,l,m,ii,jj,ef;
	PointRZ celm;

	bool fstop;

	nrec=p_nrec;
	kpnt=p_kpnt;
	krect=p_krect;
	rec=p_rec;
	pnt=p_pnt;
	rect=p_rect;
	sigma=p_sigma;
	nreg=p_nreg;
	qr=p_qr;
	qz=p_qz;
	reg=p_reg;
	rm=p_rm;
	zm=p_zm;

	if(nrec==0)return 0;

	if ((int)RecToElem.size())RecToElem.clear();
	if (OutRecList){delete[] OutRecList; OutRecList = NULL;}

	RecToElem.resize(nrec);
	OutRecList=new Output2D[nrec];

	fstop=false;
	for(k=0;k<nrec;k++)
	{
		i=FindIntervalInDoubleMas(rm,qr,rec[k].r);
		j=FindIntervalInDoubleMas(zm,qz,rec[k].z);
		if(i!=-1 && j!=-1)
		{
			m=j*(qr-1)+i;
			l=reg[m]-1;
			RecToElem[k]=l;

			celm.r=0.0;
			celm.z=0.0;
			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].oror.lem[ii]=-1;
				OutRecList[k].oz.lem[ii]=-1;
				jj=rect[l].nodes[ii]-1;
				celm.r+=pnt[jj].r;
				celm.z+=pnt[jj].z;
			}
			celm.r*=0.25;
			celm.z*=0.25;

			ef=!(rec[k].r>celm.r);
			GetSmoothForR(i,j,k,ef,l,0,rec[k],OutRecList[k].oror.lem,OutRecList[k].oror.coef,OutRecList[k].oror.v,OutRecList[k].oror.ves,OutRecList[k].oror.fl,0);

			ef=!(rec[k].z>celm.z);
			GetSmoothForZ(i,j,k,ef,l,0,rec[k],OutRecList[k].oz.lem,OutRecList[k].oz.coef,OutRecList[k].oz.v,OutRecList[k].oz.ves,OutRecList[k].oz.fl,0);
		}
		else
		{
			logfile<<"No element for reciver "<<k+1<<" r= "<<rec[k].r<<" z= "<<rec[k].z<<'\n';
			fstop=true;
		}
	}

	return fstop;
}
// Smoothing E values in receivers
void SmoothOutput2D::OutputE(double *_v2s,double *_v2c,double *results,double *resultc,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,l,m,it,nn[2];
	double s[2],ls[2],v[2],hr,hz,psi,eta;
	int ipls;
	double *v2s,*v2c;

	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceE[k/NUMBEROFABDIPOLES];
		v2s=_v2s+kpnt*ipls;
		v2c=_v2c+kpnt*ipls;

		if(OutType && OutRecList[k].oror.fl)
		{
			s[0]=0.0;
			s[1]=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oror.ves) || (it && OutRecList[k].oror.ves!=1.0))
				{					
					ls[0]=0.0;
					ls[1]=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oror.lem[j+it];
						
						v[0]=0.0;
						v[1]=0.0;
						for(l=0;l<4;l++)
						{
							v[0]+=OutRecList[k].oror.coef[j+it][l]*v2s[rect[m].nodes[l]-1];
							v[1]+=OutRecList[k].oror.coef[j+it][l]*v2c[rect[m].nodes[l]-1];
						}	
						
						ls[0]+=OutRecList[k].oror.v[it][j]*v[0];
						ls[1]+=OutRecList[k].oror.v[it][j]*v[1];
					}

					s[0]+=ls[0]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
					s[1]+=ls[1]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s[0]=(1.0-psi)*(1.0-eta)*v2s[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2s[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2s[rect[i].nodes[2]-1]+psi*eta*v2s[rect[i].nodes[3]-1];
			s[1]=(1.0-psi)*(1.0-eta)*v2c[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2c[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2c[rect[i].nodes[2]-1]+psi*eta*v2c[rect[i].nodes[3]-1];
		}

		results[k]=s[0];
		resultc[k]=s[1];
	}
}
// Smoothing B values in receivers
void SmoothOutput2D::OutputB(double *_v2s,double *_v2c,double *resultxs,double *resultxc,double *resultys,double *resultyc,
	double *resultzs,double *resultzc,vector<PointXYZ> &Ic,double *sinfi,double *cosfi,int OutType,vector<int> &RecToSourceB)
{
	double dAfidr[2],dAfidz[2],drx[2],dry[2],pr,drdx,drdy,dc=0.01;
	int i,j,k,m,it,nn[2];
	double s[2],ls[2],v[2],hr,hz,psi,eta;
	double Ix,Iy,Iz;
	int ipls;
	double *v2s,*v2c;

	for(k=0;k<nrec;k++)
	{
		j=k/NUMBEROFABDIPOLES;
		Ix=Ic[j].x;
		Iy=Ic[j].y;
		Iz=Ic[j].z;
		ipls=RecToSourceB[j];
		v2s=_v2s+kpnt*ipls;
		v2c=_v2c+kpnt*ipls;

		if(OutType && OutRecList[k].oror.fl)
		{
			s[0]=0.0;
			s[1]=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oror.ves) || (it && OutRecList[k].oror.ves!=1.0))
				{
					ls[0]=0.0;
					ls[1]=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oror.lem[j+it];

						nn[0]=rect[m].nodes[0]-1;
						nn[1]=rect[m].nodes[3]-1;

						hr=pnt[nn[1]].r-pnt[nn[0]].r;
						hz=pnt[nn[1]].z-pnt[nn[0]].z;
						eta=(rec[k].z-pnt[nn[0]].z)/hz;

						v[0]=(-(1.0-eta)*v2s[rect[m].nodes[0]-1]+(1.0-eta)*v2s[rect[m].nodes[1]-1]
							-eta*v2s[rect[m].nodes[2]-1]+eta*v2s[rect[m].nodes[3]-1])/hr;
						v[1]=(-(1.0-eta)*v2c[rect[m].nodes[0]-1]+(1.0-eta)*v2c[rect[m].nodes[1]-1]
							-eta*v2c[rect[m].nodes[2]-1]+eta*v2c[rect[m].nodes[3]-1])/hr;

						ls[0]+=OutRecList[k].oror.v[it][j]*v[0];
						ls[1]+=OutRecList[k].oror.v[it][j]*v[1];
					}
					s[0]+=ls[0]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
					s[1]+=ls[1]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s[0]=(-(1.0-eta)*v2s[rect[i].nodes[0]-1]+(1.0-eta)*v2s[rect[i].nodes[1]-1]
				-eta*v2s[rect[i].nodes[2]-1]+eta*v2s[rect[i].nodes[3]-1])/hr;
			s[1]=(-(1.0-eta)*v2c[rect[i].nodes[0]-1]+(1.0-eta)*v2c[rect[i].nodes[1]-1]
				-eta*v2c[rect[i].nodes[2]-1]+eta*v2c[rect[i].nodes[3]-1])/hr;
		}

		dAfidr[0]=s[0];
		dAfidr[1]=s[1];

		if(OutType && OutRecList[k].oz.fl)
		{
			s[0]=0.0;
			s[1]=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{
					ls[0]=0.0;
					ls[1]=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];

						nn[0]=rect[m].nodes[0]-1;
						nn[1]=rect[m].nodes[3]-1;

						hr=pnt[nn[1]].r-pnt[nn[0]].r;
						hz=pnt[nn[1]].z-pnt[nn[0]].z;
						psi=(rec[k].r-pnt[nn[0]].r)/hr;

						v[0]=(-(1.0-psi)*v2s[rect[m].nodes[0]-1]-psi*v2s[rect[m].nodes[1]-1]
								+(1.0-psi)*v2s[rect[m].nodes[2]-1]+psi*v2s[rect[m].nodes[3]-1])/hz;
						v[1]=(-(1.0-psi)*v2c[rect[m].nodes[0]-1]-psi*v2c[rect[m].nodes[1]-1]
								+(1.0-psi)*v2c[rect[m].nodes[2]-1]+psi*v2c[rect[m].nodes[3]-1])/hz;

						ls[0]+=OutRecList[k].oz.v[it][j]*v[0];
						ls[1]+=OutRecList[k].oz.v[it][j]*v[1];
					}
					s[0]+=ls[0]*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
					s[1]+=ls[1]*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;

			s[0]=(-(1.0-psi)*v2s[rect[i].nodes[0]-1]-psi*v2s[rect[i].nodes[1]-1]
				+(1.0-psi)*v2s[rect[i].nodes[2]-1]+psi*v2s[rect[i].nodes[3]-1])/hz;
			s[1]=(-(1.0-psi)*v2c[rect[i].nodes[0]-1]-psi*v2c[rect[i].nodes[1]-1]
				+(1.0-psi)*v2c[rect[i].nodes[2]-1]+psi*v2c[rect[i].nodes[3]-1])/hz;
		}

		dAfidz[0]=s[0];
		dAfidz[1]=s[1];

		pr=rec[k].r;
		if(pr<10.0*dc)pr=10.0*dc;

		drx[0]=sqrt(pr*pr+dc*dc+2*pr*dc*cosfi[k])-pr;
		drx[1]=pr-sqrt(pr*pr+dc*dc-2*pr*dc*cosfi[k]);
		dry[0]=sqrt(pr*pr+dc*dc+2*pr*dc*sinfi[k])-pr;
		dry[1]=pr-sqrt(pr*pr+dc*dc-2*pr*dc*sinfi[k]);

		drdx=(drx[0]+drx[1])/(2.0*dc);
		drdy=(dry[0]+dry[1])/(2.0*dc);

		resultxs[k]=Iz*drdy*dAfidr[0]-Iy*dAfidz[0];
		resultys[k]=Ix*dAfidz[0]-Iz*drdx*dAfidr[0];
		resultzs[k]=(Iy*drdx-Ix*drdy)*dAfidr[0];

		resultxc[k]=Iz*drdy*dAfidr[1]-Iy*dAfidz[1];
		resultyc[k]=Ix*dAfidz[1]-Iz*drdx*dAfidr[1];
		resultzc[k]=(Iy*drdx-Ix*drdy)*dAfidr[1];
	}
}

void SmoothOutput2D::GetSmoothForR(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int rlem[4],double rcoef[4][4],double rv[2][3],double &ves,int &fl,int stype)
{
	int n,ii,kk,ll,lt,mm,tt;
	double cr[3],r0,r3,z0,z3,invhrhz;
	bool fpp[2];
	int lef;
	PointRZ clem,lrecz;

	fl=0;
	ves=0.0;

	for(kk=0;kk<4;kk++)rlem[kk]=-1;

	rlem[1+ef]=l;

	if(stype)
	{
		clem.r=0.0;
		clem.z=0.0;
		for(kk=0;kk<4;kk++)
		{
			clem.r+=pnt[rect[l].nodes[kk]-1].r;
			clem.z+=pnt[rect[l].nodes[kk]-1].z;
		}
		clem.r*=0.25;
		clem.z*=0.25;

		lrecz.r=clem.r;
		lrecz.z=reck.z;
		lef=!(lrecz.z>clem.z);

		tt=1+ef;
		OutRecList[k].llemz[tt]=1+lef;
		OutRecList[k].lrecz[tt]=lrecz;
		GetSmoothForZ(i,j,k,lef,l,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
	}

	lt=l;
	n=0;
	for(ii=i-1;ii>=0 && n<(1+ef);ii--)
	{
		mm=j*(qr-1)+ii;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			rlem[(1+ef)-n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecz.r=clem.r;
				lrecz.z=reck.z;
				lef=!(lrecz.z>clem.z);

				tt=(1+ef)-n;
				OutRecList[k].llemz[tt]=1+lef;
				OutRecList[k].lrecz[tt]=lrecz;
				GetSmoothForZ(ii,j,k,lef,ll,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
			}
		}
	}

	if(ef && rlem[1]==-1)
	{
		rlem[1]=rlem[1+ef];
		rlem[1+ef]=-1;

		if(stype)
		{
			OutRecList[k].llemz[1]=OutRecList[k].llemz[1+ef];
			OutRecList[k].lrecz[1]=OutRecList[k].lrecz[1+ef];
			OutRecList[k].ozfr[1].fl=OutRecList[k].ozfr[1+ef].fl;
			OutRecList[k].ozfr[1].ves=OutRecList[k].ozfr[1+ef].ves;

			for(ii=0;ii<2;ii++)
			{
				for(kk=0;kk<3;kk++)
				{
					OutRecList[k].ozfr[1].v[ii][kk]=OutRecList[k].ozfr[1+ef].v[ii][kk];
				}
			}

			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].ozfr[1].lem[ii]=OutRecList[k].ozfr[1+ef].lem[ii];
				if(OutRecList[k].ozfr[1].lem[ii]!=-1)
				{
					for(kk=0;kk<4;kk++)
					{
						OutRecList[k].ozfr[1].coef[ii][kk]=OutRecList[k].ozfr[1+ef].coef[ii][kk];
					}
				}
			}
		}
		ef=0;
	}

	lt=l;
	n=0;
	for(ii=i+1;ii<qr-1 && n<(2-ef);ii++)
	{
		mm=j*(qr-1)+ii;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			rlem[(1+ef)+n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecz.r=clem.r;
				lrecz.z=reck.z;
				lef=!(lrecz.z>clem.z);

				tt=(1+ef)+n;
				OutRecList[k].llemz[tt]=1+lef;
				OutRecList[k].lrecz[tt]=lrecz;
				GetSmoothForZ(ii,j,k,lef,ll,0,lrecz,OutRecList[k].ozfr[tt].lem,OutRecList[k].ozfr[tt].coef,OutRecList[k].ozfr[tt].v,OutRecList[k].ozfr[tt].ves,OutRecList[k].ozfr[tt].fl,0);
			}
		}
	}

	fpp[0]=(rlem[0]!=-1 && rlem[1]!=-1 && rlem[2]!=-1);
	fpp[1]=(rlem[1]!=-1 && rlem[2]!=-1 && rlem[3]!=-1);

	if(fpp[0])
	{
		for(kk=0;kk<3;kk++)
		{
			cr[kk]=0.0;
			ll=rlem[kk];
			for(ii=0;ii<4;ii++)
			{
				cr[kk]+=pnt[rect[ll].nodes[ii]-1].r;
			}
			cr[kk]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			rcoef[kk][0]=((r3-cr[kk])*(z3-reck.z))*invhrhz;
			rcoef[kk][1]=((cr[kk]-r0)*(z3-reck.z))*invhrhz;
			rcoef[kk][2]=((r3-cr[kk])*(reck.z-z0))*invhrhz;
			rcoef[kk][3]=((cr[kk]-r0)*(reck.z-z0))*invhrhz;
		}
		if(!fpp[1])
		{
			fl=1;
			ves=1.0;
		}
		else
		{
			fl=2;
			ves=1.0-(reck.r-cr[1])/(cr[2]-cr[1]);
		}
		rv[0][0]=((reck.r-cr[1])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[2]-cr[0]));
		rv[0][1]=((reck.r-cr[0])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[1]-cr[2]));
		rv[0][2]=((reck.r-cr[0])*(reck.r-cr[1]))/((cr[2]-cr[0])*(cr[2]-cr[1]));
	}

	if(fpp[1])
	{
		if(!fpp[0])
		{
			for(kk=0;kk<3;kk++)
			{
				cr[kk]=0.0;
				ll=rlem[kk+1];
				for(ii=0;ii<4;ii++)
				{
					cr[kk]+=pnt[rect[ll].nodes[ii]-1].r;
				}
				cr[kk]*=0.25;
				
				r0=pnt[rect[ll].nodes[0]-1].r;
				z0=pnt[rect[ll].nodes[0]-1].z;
				r3=pnt[rect[ll].nodes[3]-1].r;
				z3=pnt[rect[ll].nodes[3]-1].z;

				invhrhz=1.0/((r3-r0)*(z3-z0));

				rcoef[kk+1][0]=((r3-cr[kk])*(z3-reck.z))*invhrhz;
				rcoef[kk+1][1]=((cr[kk]-r0)*(z3-reck.z))*invhrhz;
				rcoef[kk+1][2]=((r3-cr[kk])*(reck.z-z0))*invhrhz;
				rcoef[kk+1][3]=((cr[kk]-r0)*(reck.z-z0))*invhrhz;
			}
			fl=1;
			ves=0.0;
		}
		else
		{
			cr[0]=cr[1];
			cr[1]=cr[2];
			cr[2]=0.0;
			ll=rlem[3];
			for(ii=0;ii<4;ii++)
			{
				cr[2]+=pnt[rect[ll].nodes[ii]-1].r;
			}
			cr[2]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			rcoef[3][0]=((r3-cr[2])*(z3-reck.z))*invhrhz;
			rcoef[3][1]=((cr[2]-r0)*(z3-reck.z))*invhrhz;
			rcoef[3][2]=((r3-cr[2])*(reck.z-z0))*invhrhz;
			rcoef[3][3]=((cr[2]-r0)*(reck.z-z0))*invhrhz;
		}
		rv[1][0]=((reck.r-cr[1])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[2]-cr[0]));
		rv[1][1]=((reck.r-cr[0])*(reck.r-cr[2]))/((cr[1]-cr[0])*(cr[1]-cr[2]));
		rv[1][2]=((reck.r-cr[0])*(reck.r-cr[1]))/((cr[2]-cr[0])*(cr[2]-cr[1]));
	}

	if(stype)
	{
		for(kk=0;kk<4;kk++)
		{
			if(rlem[kk]==-1)OutRecList[k].ozfr[kk].fl=0;
		}
	}
}

void SmoothOutput2D::GetSmoothForZ(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int zlem[4],double zcoef[4][4],double zv[2][3],double &ves,int &fl,int stype)
{
	int n,ii,kk,ll,lt,mm,tt;
	double cz[3],r0,r3,z0,z3,invhrhz;
	bool fpp[2];
	int lef;
	PointRZ clem,lrecr;

	fl=0;
	ves=0.0;

	for(kk=0;kk<4;kk++)zlem[kk]=-1;

	zlem[1+ef]=l;

	if(stype)
	{
		clem.r=0.0;
		clem.z=0.0;
		for(kk=0;kk<4;kk++)
		{
			clem.r+=pnt[rect[l].nodes[kk]-1].r;
			clem.z+=pnt[rect[l].nodes[kk]-1].z;
		}
		clem.r*=0.25;
		clem.z*=0.25;

		lrecr.r=reck.r;
		lrecr.z=clem.z;
		lef=!(lrecr.r>clem.r);

		tt=1+ef;
		OutRecList[k].llemr[tt]=1+lef;
		OutRecList[k].lrecr[tt]=lrecr;
		GetSmoothForR(i,j,k,lef,l,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
	}

	lt=l;
	n=0;
	for(ii=j-1;ii>=0 && n<(1+ef);ii--)
	{
		mm=ii*(qr-1)+i;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			zlem[(1+ef)-n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecr.r=reck.r;
				lrecr.z=clem.z;
				lef=!(lrecr.r>clem.r);

				tt=(1+ef)-n;
				OutRecList[k].llemr[tt]=1+lef;
				OutRecList[k].lrecr[tt]=lrecr;
				GetSmoothForR(i,ii,k,lef,ll,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
			}
		}
	}

	if(ef && zlem[1]==-1)
	{
		zlem[1]=zlem[1+ef];
		zlem[1+ef]=-1;
		
		if(stype)
		{
			OutRecList[k].llemr[1]=OutRecList[k].llemr[1+ef];
			OutRecList[k].lrecr[1]=OutRecList[k].lrecr[1+ef];
			OutRecList[k].orfz[1].fl=OutRecList[k].orfz[1+ef].fl;
			OutRecList[k].orfz[1].ves=OutRecList[k].orfz[1+ef].ves;

			for(ii=0;ii<2;ii++)
			{
				for(kk=0;kk<3;kk++)
				{
					OutRecList[k].orfz[1].v[ii][kk]=OutRecList[k].orfz[1+ef].v[ii][kk];
				}
			}

			for(ii=0;ii<4;ii++)
			{
				OutRecList[k].orfz[1].lem[ii]=OutRecList[k].orfz[1+ef].lem[ii];
				if(OutRecList[k].orfz[1].lem[ii]!=-1)
				{
					for(kk=0;kk<4;kk++)
					{
						OutRecList[k].orfz[1].coef[ii][kk]=OutRecList[k].orfz[1+ef].coef[ii][kk];
					}
				}
			}
		}
		ef=0;
	}

	lt=l;
	n=0;
	for(ii=j+1;ii<qz-1 && n<(2-ef);ii++)
	{
		mm=ii*(qr-1)+i;
		ll=reg[mm]-1;
		if(ll!=lt)
		{
			if(MatFlag && rect[l].mtr!=rect[ll].mtr)break;
			n+=1;
			zlem[(1+ef)+n]=ll;
			lt=ll;

			if(stype)
			{
				clem.r=0.0;
				clem.z=0.0;
				for(kk=0;kk<4;kk++)
				{
					clem.r+=pnt[rect[ll].nodes[kk]-1].r;
					clem.z+=pnt[rect[ll].nodes[kk]-1].z;
				}
				clem.r*=0.25;
				clem.z*=0.25;

				lrecr.r=reck.r;
				lrecr.z=clem.z;
				lef=!(lrecr.r>clem.r);

				tt=(1+ef)+n;
				OutRecList[k].llemr[tt]=1+lef;
				OutRecList[k].lrecr[tt]=lrecr;
				GetSmoothForR(i,ii,k,lef,ll,0,lrecr,OutRecList[k].orfz[tt].lem,OutRecList[k].orfz[tt].coef,OutRecList[k].orfz[tt].v,OutRecList[k].orfz[tt].ves,OutRecList[k].orfz[tt].fl,0);
			}
		}
	}

	fpp[0]=(zlem[0]!=-1 && zlem[1]!=-1 && zlem[2]!=-1);
	fpp[1]=(zlem[1]!=-1 && zlem[2]!=-1 && zlem[3]!=-1);

	if(fpp[0])
	{
		for(kk=0;kk<3;kk++)
		{
			cz[kk]=0.0;
			ll=zlem[kk];
			for(ii=0;ii<4;ii++)
			{
				cz[kk]+=pnt[rect[ll].nodes[ii]-1].z;
			}
			cz[kk]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			zcoef[kk][0]=((r3-reck.r)*(z3-cz[kk]))*invhrhz;
			zcoef[kk][1]=((reck.r-r0)*(z3-cz[kk]))*invhrhz;
			zcoef[kk][2]=((r3-reck.r)*(cz[kk]-z0))*invhrhz;
			zcoef[kk][3]=((reck.r-r0)*(cz[kk]-z0))*invhrhz;
		}
		if(!fpp[1])
		{
			fl=1;
			ves=1.0;
		}
		else
		{
			fl=2;
			ves=1.0-(reck.z-cz[1])/(cz[2]-cz[1]);
		}
		zv[0][0]=((reck.z-cz[1])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[2]-cz[0]));
		zv[0][1]=((reck.z-cz[0])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[1]-cz[2]));
		zv[0][2]=((reck.z-cz[0])*(reck.z-cz[1]))/((cz[2]-cz[0])*(cz[2]-cz[1]));
	}

	if(fpp[1])
	{
		if(!fpp[0])
		{
			for(kk=0;kk<3;kk++)
			{
				cz[kk]=0.0;
				ll=zlem[kk+1];
				for(ii=0;ii<4;ii++)
				{
					cz[kk]+=pnt[rect[ll].nodes[ii]-1].z;
				}
				cz[kk]*=0.25;
				
				r0=pnt[rect[ll].nodes[0]-1].r;
				z0=pnt[rect[ll].nodes[0]-1].z;
				r3=pnt[rect[ll].nodes[3]-1].r;
				z3=pnt[rect[ll].nodes[3]-1].z;

				invhrhz=1.0/((r3-r0)*(z3-z0));

				zcoef[kk+1][0]=((r3-reck.r)*(z3-cz[kk]))*invhrhz;
				zcoef[kk+1][1]=((reck.r-r0)*(z3-cz[kk]))*invhrhz;
				zcoef[kk+1][2]=((r3-reck.r)*(cz[kk]-z0))*invhrhz;
				zcoef[kk+1][3]=((reck.r-r0)*(cz[kk]-z0))*invhrhz;
			}
			fl=1;
			ves=0.0;
		}
		else
		{
			cz[0]=cz[1];
			cz[1]=cz[2];
			cz[2]=0.0;
			ll=zlem[3];
			for(ii=0;ii<4;ii++)
			{
				cz[2]+=pnt[rect[ll].nodes[ii]-1].z;
			}
			cz[2]*=0.25;
			
			r0=pnt[rect[ll].nodes[0]-1].r;
			z0=pnt[rect[ll].nodes[0]-1].z;
			r3=pnt[rect[ll].nodes[3]-1].r;
			z3=pnt[rect[ll].nodes[3]-1].z;

			invhrhz=1.0/((r3-r0)*(z3-z0));

			zcoef[3][0]=((r3-reck.r)*(z3-cz[2]))*invhrhz;
			zcoef[3][1]=((reck.r-r0)*(z3-cz[2]))*invhrhz;
			zcoef[3][2]=((r3-reck.r)*(cz[2]-z0))*invhrhz;
			zcoef[3][3]=((reck.r-r0)*(cz[2]-z0))*invhrhz;
		}
		zv[1][0]=((reck.z-cz[1])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[2]-cz[0]));
		zv[1][1]=((reck.z-cz[0])*(reck.z-cz[2]))/((cz[1]-cz[0])*(cz[1]-cz[2]));
		zv[1][2]=((reck.z-cz[0])*(reck.z-cz[1]))/((cz[2]-cz[0])*(cz[2]-cz[1]));
	}

	if(stype)
	{
		for(kk=0;kk<4;kk++)
		{
			if(zlem[kk]==-1)OutRecList[k].orfz[kk].fl=0;
		}
	}
}

void SmoothOutput2D::OutputE(double *v2s,double *v2c,double *results,double *resultc,int OutType)
{
	int i,j,k,l,m,it,nn[2];
	double s[2],ls[2],v[2],hr,hz,psi,eta;

	for(k=0;k<nrec;k++)
	{
		if(OutType && OutRecList[k].oror.fl)
		{
			s[0]=0.0;
			s[1]=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oror.ves) || (it && OutRecList[k].oror.ves!=1.0))
				{					
					ls[0]=0.0;
					ls[1]=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oror.lem[j+it];
						
						v[0]=0.0;
						v[1]=0.0;
						for(l=0;l<4;l++)
						{
							v[0]+=OutRecList[k].oror.coef[j+it][l]*v2s[rect[m].nodes[l]-1];
							v[1]+=OutRecList[k].oror.coef[j+it][l]*v2c[rect[m].nodes[l]-1];
						}	
						
						ls[0]+=OutRecList[k].oror.v[it][j]*v[0];
						ls[1]+=OutRecList[k].oror.v[it][j]*v[1];
					}

					s[0]+=ls[0]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
					s[1]+=ls[1]*((it)? 1.0-OutRecList[k].oror.ves : OutRecList[k].oror.ves);
				}
			}
		}
		else
		{
			i=RecToElem[k];

			nn[0]=rect[i].nodes[0]-1;
			nn[1]=rect[i].nodes[3]-1;

			hr=pnt[nn[1]].r-pnt[nn[0]].r;
			hz=pnt[nn[1]].z-pnt[nn[0]].z;
			psi=(rec[k].r-pnt[nn[0]].r)/hr;
			eta=(rec[k].z-pnt[nn[0]].z)/hz;

			s[0]=(1.0-psi)*(1.0-eta)*v2s[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2s[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2s[rect[i].nodes[2]-1]+psi*eta*v2s[rect[i].nodes[3]-1];
			s[1]=(1.0-psi)*(1.0-eta)*v2c[rect[i].nodes[0]-1]+psi*(1.0-eta)*v2c[rect[i].nodes[1]-1]
				+(1.0-psi)*eta*v2c[rect[i].nodes[2]-1]+psi*eta*v2c[rect[i].nodes[3]-1];
		}

		results[k]=s[0];
		resultc[k]=s[1];
	}
}
