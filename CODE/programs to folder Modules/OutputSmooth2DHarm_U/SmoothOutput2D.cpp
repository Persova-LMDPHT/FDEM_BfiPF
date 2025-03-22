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
 *  This file contains code for direct FEM output E and B fields in receivers.
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
extern bool UseRotU;
extern bool UseGradDivA;

SmoothOutput2D::SmoothOutput2D()
{
	nrec=kpnt=krect=0;
	rec=NULL;
	pnt=NULL;
	rect=NULL;
	OutRecList=NULL;
}

SmoothOutput2D::~SmoothOutput2D()
{
	if(OutRecList){delete [] OutRecList; OutRecList=NULL;}
}

// Initialization of smoothing output (find elements for receivers)
int SmoothOutput2D::Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
						 double *p_sigmaZ,double *p_dpr,int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm,double _w)
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
	sigmaZ=p_sigmaZ;
	dpr=p_dpr;
	nreg=p_nreg;
	qr=p_qr;
	qz=p_qz;
	reg=p_reg;
	rm=p_rm;
	zm=p_zm;
	w=_w;

	if(nrec==0)return 0;

	RecToElem.resize(nrec);
	OutRecList=new Output2D[nrec];

	fstop=false;
	#pragma omp parallel reduction(||: fstop)
	{
		#pragma omp for private(i,j,k,l,m,ii,jj,ef,celm)
		for (k = 0; k < nrec; k++)
		{
			i = FindIntervalInDoubleMas(rm, qr, rec[k].r);
			j = FindIntervalInDoubleMas(zm, qz, rec[k].z);
			if (i != -1 && j != -1)
			{
				m = j * (qr - 1) + i;
				l = reg[m] - 1;

				if (matrec)
				{
					if (matrec[k] && rect[l].mtr != matrec[k])
					{
						jj = 0;
						if (fabs(rec[k].z - pnt[rect[l].nodes[0] - 1].z) < 1e-6)
						{
							if (j > 0)
							{
								m = (j - 1) * (qr - 1) + i;
								l = reg[m] - 1;
								if (rect[l].mtr != matrec[k])
								{
									jj = 1;
								}
								else
								{
									j--;
								}
							}
							else
							{
								jj = 1;
							}
						}
						else if (fabs(rec[k].z - pnt[rect[l].nodes[3] - 1].z) < 1e-6)
						{
							if (j < qz - 1)
							{
								m = (j + 1) * (qr - 1) + i;
								l = reg[m] - 1;
								if (rect[l].mtr != matrec[k])
								{
									jj = 1;
								}
								else
								{
									j++;
								}
							}
							else
							{
								jj = 1;
							}
						}
						else
						{
							jj = 1;
						}

						if (jj)
						{
							logfile << "Material problem for reciver " << k + 1 << " r= " << rec[k].r << " z= " << rec[k].z << '\n';
							fstop = true;
						}
					}
				}

				if(!fstop)
				{
					RecToElem[k] = l;

					celm.r = 0.0;
					celm.z = 0.0;
					for (ii = 0; ii < 4; ii++)
					{
						OutRecList[k].orobj.lem[ii] = -1;
						OutRecList[k].oz.lem[ii] = -1;
						jj = rect[l].nodes[ii] - 1;
						celm.r += pnt[jj].r;
						celm.z += pnt[jj].z;
					}
					celm.r *= 0.25;
					celm.z *= 0.25;

					// r
					ef = !(rec[k].r > celm.r);
					GetSmoothForR(i, j, k, ef, l, 0, rec[k], OutRecList[k].orobj.lem, OutRecList[k].orobj.coef, OutRecList[k].orobj.v, OutRecList[k].orobj.ves, OutRecList[k].orobj.fl, 1);

					// z
					ef = !(rec[k].z > celm.z);
					GetSmoothForZ(i, j, k, ef, l, 0, rec[k], OutRecList[k].oz.lem, OutRecList[k].oz.coef, OutRecList[k].oz.v, OutRecList[k].oz.ves, OutRecList[k].oz.fl, 1);
				}
			}
			else
			{
				logfile << "No element for reciver " << k + 1 << " r= " << rec[k].r << " z= " << rec[k].z << '\n';
				fstop = true;
			}
		}
	}

	return fstop;
}
// Smoothing Er values in receivers
void SmoothOutput2D::Output_U_Er(double *_v2us,double *_v2uc,double *_v2axs,double *_v2axc,double *results,double *resultc,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,m,it,node0,node1,node2,node3,ipls;
	double ss,lss,vs,sc,lsc,vc,hr,hz,psi,eta,znam,sig,weps,fs,fc;
	double *v2us,*v2uc,*v2axs,*v2axc;

	#pragma omp parallel for private(i,j,k,m,it,node0,node1,node2,node3,ss,sc,lss,lsc,vs,vc,hr,hz,psi,eta,znam,sig,weps,fs,fc,ipls,v2us,v2uc,v2axs,v2axc)
	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceE[k/2];
		v2us=_v2us+kpnt*ipls;
		v2uc=_v2uc+kpnt*ipls;
		v2axs=_v2axs+kpnt*ipls;
		v2axc=_v2axc+kpnt*ipls;

		if(OutType && OutRecList[k].orobj.fl)
		{
			ss=sc=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].orobj.ves) || (it && OutRecList[k].orobj.ves!=1.0))
				{
					lss=lsc=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].orobj.lem[j+it];

						node0=rect[m].nodes[0]-1;
						node1=rect[m].nodes[1]-1;
						node2=rect[m].nodes[2]-1;
						node3=rect[m].nodes[3]-1;

						hr=pnt[node3].r-pnt[node0].r;
						hz=pnt[node3].z-pnt[node0].z;
				
						GetSmoothVarForR_U(k,j+it,v2us,v2uc,v2axs,v2axc,vs,vc,OutType);
						lss+=OutRecList[k].orobj.v[it][j]*vs;
						lsc+=OutRecList[k].orobj.v[it][j]*vc;
					}

					ss+=lss*((it)? 1.0-OutRecList[k].orobj.ves : OutRecList[k].orobj.ves);
					sc+=lsc*((it)? 1.0-OutRecList[k].orobj.ves : OutRecList[k].orobj.ves);
				}
			}
		}
		else
		{
			ss=sc=0.0;
			i=RecToElem[k];

			node0=rect[i].nodes[0]-1;
			node1=rect[i].nodes[1]-1;
			node2=rect[i].nodes[2]-1;
			node3=rect[i].nodes[3]-1;

			hr=pnt[node3].r-pnt[node0].r;
			hz=pnt[node3].z-pnt[node0].z;
			
			psi=(rec[k].r-pnt[node0].r)/hr;
			eta=(rec[k].z-pnt[node0].z)/hz;

			if(UseRotU)
			{
				ss-=(-(1.0-psi)*v2us[node0]-psi*v2us[node1]+(1.0-psi)*v2us[node2]+psi*v2us[node3])/hz;
				sc-=(-(1.0-psi)*v2uc[node0]-psi*v2uc[node1]+(1.0-psi)*v2uc[node2]+psi*v2uc[node3])/hz;
			}

			if(UseGradDivA)
			{
				ss+=(-(1.0-eta)*v2axs[node0]+(1.0-eta)*v2axs[node1]-eta*v2axs[node2]+eta*v2axs[node3])/hr;
				sc+=(-(1.0-eta)*v2axc[node0]+(1.0-eta)*v2axc[node1]-eta*v2axc[node2]+eta*v2axc[node3])/hr;
			}

			sig=sigma[rect[i].mtr-1];
			weps=w*dpr[rect[i].mtr-1];

			znam=sig*sig+weps*weps;

			fs=ss;
			fc=sc;

			ss=(sig*fs+weps*fc)/znam;
			sc=(-weps*fs+sig*fc)/znam;
		}

		results[k]=ss/MU0;
		resultc[k]=sc/MU0;
	}
}
// Smoothing Ez values in receivers
void SmoothOutput2D::Output_U_Ez(double *_v2us,double *_v2uc,double *_v2axs,double *_v2axc,double *results,double *resultc,int OutType,vector<int> &RecToSourceE)
{
	int i,j,k,m,it,node0,node1,node2,node3,ipls;
	double ss,lss,vs,sc,lsc,vc,hr,hz,psi,eta,r;
	double *v2us,*v2uc,*v2axs,*v2axc;
	double sig,weps,znam;
	
	#pragma omp parallel for private(i,j,k,m,it,node0,node1,node2,node3,ss,sc,lss,lsc,vs,vc,hr,hz,psi,eta,r,sig,weps,znam,ipls,v2us,v2uc,v2axs,v2axc)
	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceE[k/2];
		v2us=_v2us+kpnt*ipls;
		v2uc=_v2uc+kpnt*ipls;
		v2axs=_v2axs+kpnt*ipls;
		v2axc=_v2axc+kpnt*ipls;

		if(OutType && OutRecList[k].oz.fl)
		{
			ss=sc=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{					
					lss=lsc=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];

						vs=GetSmoothValForZ_U(k,j+it,v2us,v2axs,OutType);
						lss+=OutRecList[k].oz.v[it][j]*vs;
						vc=GetSmoothValForZ_U(k,j+it,v2uc,v2axc,OutType);
						lsc+=OutRecList[k].oz.v[it][j]*vc;
					}

					ss+=lss*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
					sc+=lsc*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			ss=sc=0.0;

			if (OutType && OutRecList[k].orobj.fl)
			{
				ss = sc = 0.0;
				for (it = 0; it < 2; it++)
				{
					if ((!it && OutRecList[k].orobj.ves) || (it && OutRecList[k].orobj.ves != 1.0))
					{
						lss = lsc = 0.0;
						for (j = 0; j < 3; j++)
						{
							m = OutRecList[k].orobj.lem[j + it];

							node0 = rect[m].nodes[0] - 1;
							node1 = rect[m].nodes[1] - 1;
							node2 = rect[m].nodes[2] - 1;
							node3 = rect[m].nodes[3] - 1;

							hr = pnt[node3].r - pnt[node0].r;
							hz = pnt[node3].z - pnt[node0].z;

							GetSmoothVarForZR_U(k, j + it, v2us, v2uc, v2axs, v2axc, vs, vc);
							lss += OutRecList[k].orobj.v[it][j] * vs;
							lsc += OutRecList[k].orobj.v[it][j] * vc;
						}

						ss += lss * ((it) ? 1.0 - OutRecList[k].orobj.ves : OutRecList[k].orobj.ves);
						sc += lsc * ((it) ? 1.0 - OutRecList[k].orobj.ves : OutRecList[k].orobj.ves);
					}
				}
			}
			else
			{
				i=RecToElem[k];

				node0=rect[i].nodes[0]-1;
				node1=rect[i].nodes[1]-1;
				node2=rect[i].nodes[2]-1;
				node3=rect[i].nodes[3]-1;

				hr=pnt[node3].r-pnt[node0].r;
				hz=pnt[node3].z-pnt[node0].z;
			
				psi=(rec[k].r-pnt[node0].r)/hr;
				eta=(rec[k].z-pnt[node0].z)/hz;

				r=rec[k].r;

				if(UseRotU)
				{
					ss+=((1.0-psi)*(1.0-eta)*v2us[node0]+psi*(1.0-eta)*v2us[node1]+(1.0-psi)*eta*v2us[node2]+psi*eta*v2us[node3])/r+
						(-(1.0-eta)*v2us[node0]+(1.0-eta)*v2us[node1]-eta*v2us[node2]+eta*v2us[node3])/hr;
					sc+=((1.0-psi)*(1.0-eta)*v2uc[node0]+psi*(1.0-eta)*v2uc[node1]+(1.0-psi)*eta*v2uc[node2]+psi*eta*v2uc[node3])/r+
						(-(1.0-eta)*v2uc[node0]+(1.0-eta)*v2uc[node1]-eta*v2uc[node2]+eta*v2uc[node3])/hr;
				}

				if(UseGradDivA)
				{
					ss+=(-(1.0-psi)*v2axs[node0]-psi*v2axs[node1]+(1.0-psi)*v2axs[node2]+psi*v2axs[node3])/hz;
					sc+=(-(1.0-psi)*v2axc[node0]-psi*v2axc[node1]+(1.0-psi)*v2axc[node2]+psi*v2axc[node3])/hz;
				}
			}
		}

		i=RecToElem[k];

		sig=sigmaZ[rect[i].mtr-1];
		weps=w*dpr[rect[i].mtr-1];

		znam=sig*sig+weps*weps;
		results[k]=(sig*ss+weps*sc)/znam/MU0;
		resultc[k]=(-weps*ss+sig*sc)/znam/MU0;
	}
}
// Smoothing Bfi values in receivers
void SmoothOutput2D::Output_U_B(double *_v2us,double *_v2uc,double *results,double *resultc,int OutType,vector<int> &RecToSourceB)
{
	int i,j,k,m,it,node0,node1,node2,node3,ipls;
	double ss,lss,vs,sc,lsc,vc,hr,hz,psi,eta,r;
	double *v2us,*v2uc;
	
	for(k=0;k<nrec;k++)
	{
		ipls=RecToSourceB[k/2];
		v2us=_v2us+kpnt*ipls;
		v2uc=_v2uc+kpnt*ipls;

		if(OutType && OutRecList[k].oz.fl)
		{
			ss=sc=0.0;
			for(it=0;it<2;it++)
			{
				if((!it && OutRecList[k].oz.ves) || (it && OutRecList[k].oz.ves!=1.0))
				{
					lss=lsc=0.0;
					for(j=0;j<3;j++)
					{
						m=OutRecList[k].oz.lem[j+it];

						vs=GetSmoothValForB_U(k,j+it,v2us,OutType);
						lss+=OutRecList[k].oz.v[it][j]*vs;
						vc=GetSmoothValForB_U(k,j+it,v2uc,OutType);
						lsc+=OutRecList[k].oz.v[it][j]*vc;
					}

					ss+=lss*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
					sc+=lsc*((it)? 1.0-OutRecList[k].oz.ves : OutRecList[k].oz.ves);
				}
			}
		}
		else
		{
			ss=sc=0.0;
			i=RecToElem[k];

			node0=rect[i].nodes[0]-1;
			node1=rect[i].nodes[1]-1;
			node2=rect[i].nodes[2]-1;
			node3=rect[i].nodes[3]-1;

			hr=pnt[node3].r-pnt[node0].r;
			hz=pnt[node3].z-pnt[node0].z;
			
			psi=(rec[k].r-pnt[node0].r)/hr;
			eta=(rec[k].z-pnt[node0].z)/hz;

			r=rec[k].r;

			ss+=((1.0-psi)*(1.0-eta)*v2us[node0]+psi*(1.0-eta)*v2us[node1]+(1.0-psi)*eta*v2us[node2]+psi*eta*v2us[node3]);
			sc+=((1.0-psi)*(1.0-eta)*v2uc[node0]+psi*(1.0-eta)*v2uc[node1]+(1.0-psi)*eta*v2uc[node2]+psi*eta*v2uc[node3]);
		}

		results[k]=ss;
		resultc[k]=sc;
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
			if(rect[ll].mtr==1 || (MatFlag && rect[l].mtr!=rect[ll].mtr))break;
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
			if(rect[ll].mtr==1 || (MatFlag && rect[l].mtr!=rect[ll].mtr))break;
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

double SmoothOutput2D::GetSmoothValForZ_U(int k,int jj,double *v2u,double *v2a,int OutType)
{
	int i,j,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta,r;

	if(OutType || !OutRecList[k].orfz[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].orfz[jj].ves) || (it && OutRecList[k].orfz[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].orfz[jj].lem[j+it];

					node0=rect[m].nodes[0]-1;
					node1=rect[m].nodes[1]-1;
					node2=rect[m].nodes[2]-1;
					node3=rect[m].nodes[3]-1;

					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					
					psi=eta=0.5;
					r=0.5*(pnt[node3].r+pnt[node0].r);				

					v=0.0;
					
					if(UseRotU)
					{
						v+=((1.0-psi)*(1.0-eta)*v2u[node0]+psi*(1.0-eta)*v2u[node1]+(1.0-psi)*eta*v2u[node2]+psi*eta*v2u[node3])/r+
							(-(1.0-eta)*v2u[node0]+(1.0-eta)*v2u[node1]-eta*v2u[node2]+eta*v2u[node3])/hr;
					}

					if(UseGradDivA)
					{
						v+=(-(1.0-psi)*v2a[node0]-psi*v2a[node1]+(1.0-psi)*v2a[node2]+psi*v2a[node3])/hz;
					}

					ls+=OutRecList[k].orfz[jj].v[it][j]*v;
				}
				s+=ls*((it)? 1.0-OutRecList[k].orfz[jj].ves : OutRecList[k].orfz[jj].ves);
			}
		}
	}
	else
	{
		s=0.0;
		i=OutRecList[k].llemr[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;

		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;
		
		psi=(rec[k].r-pnt[node0].r)/hr;
		eta=0.5;
		r=rec[k].r;

		if(UseRotU)
		{
			s+=((1.0-psi)*(1.0-eta)*v2u[node0]+psi*(1.0-eta)*v2u[node1]+(1.0-psi)*eta*v2u[node2]+psi*eta*v2u[node3])/r+
				(-(1.0-eta)*v2u[node0]+(1.0-eta)*v2u[node1]-eta*v2u[node2]+eta*v2u[node3])/hr;
		}

		if(UseGradDivA)
		{
			s+=(-(1.0-psi)*v2a[node0]-psi*v2a[node1]+(1.0-psi)*v2a[node2]+psi*v2a[node3])/hz;
		}
	}

	return s;
}

void SmoothOutput2D::GetSmoothVarForR_U(int k,int jj,double *v2us,double *v2uc,double *v2as,double *v2ac,double &rs,double &rc,int OutType)
{
	int i,j,m,it,node0,node1,node2,node3;
	double ls[2],v[2],hr,hz,psi,eta,sig,weps,znam;

	if(OutType && OutRecList[k].ozfr[jj].fl)
	{
		rs=rc=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].ozfr[jj].ves) || (it && OutRecList[k].ozfr[jj].ves!=1.0))
			{
				ls[0]=ls[1]=0.0;
				for(j=0;j<3;j++)
				{
					v[0]=v[1]=0.0;

					m=OutRecList[k].ozfr[jj].lem[j+it];

					node0=rect[m].nodes[0]-1;
					node1=rect[m].nodes[1]-1;
					node2=rect[m].nodes[2]-1;
					node3=rect[m].nodes[3]-1;
					
					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;
					
					psi=eta=0.5;

					if(UseRotU)
					{
						v[0]-=(-(1.0-psi)*v2us[node0]-psi*v2us[node1]+(1.0-psi)*v2us[node2]+psi*v2us[node3])/hz;
						v[1]-=(-(1.0-psi)*v2uc[node0]-psi*v2uc[node1]+(1.0-psi)*v2uc[node2]+psi*v2uc[node3])/hz;
					}

					if(UseGradDivA)
					{
						v[0]+=(-(1.0-eta)*v2as[node0]+(1.0-eta)*v2as[node1]-eta*v2as[node2]+eta*v2as[node3])/hr;
						v[1]+=(-(1.0-eta)*v2ac[node0]+(1.0-eta)*v2ac[node1]-eta*v2ac[node2]+eta*v2ac[node3])/hr;
					}

					sig=sigma[rect[m].mtr-1];
					weps=w*dpr[rect[m].mtr-1];

					znam=sig*sig+weps*weps;

					ls[0]+=(OutRecList[k].ozfr[jj].v[it][j])*(sig*v[0]+weps*v[1])/znam;;
					ls[1]+=(OutRecList[k].ozfr[jj].v[it][j])*(-weps*v[0]+sig*v[1])/znam;
				}
				rs+=ls[0]*((it)? 1.0-OutRecList[k].ozfr[jj].ves : OutRecList[k].ozfr[jj].ves);
				rc+=ls[1]*((it)? 1.0-OutRecList[k].ozfr[jj].ves : OutRecList[k].ozfr[jj].ves);
			}
		}
	}
	else
	{
		i=OutRecList[k].orobj.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;
		
		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;

		psi=0.5;
		eta=(rec[k].z-pnt[node0].z)/hz;

		v[0]=v[1]=0.0;

		if(UseRotU)
		{
			v[0] -= (-(1.0 - psi) * v2us[node0] - psi * v2us[node1] + (1.0 - psi) * v2us[node2] + psi * v2us[node3]) / hz;
			v[1] -= (-(1.0 - psi) * v2uc[node0] - psi * v2uc[node1] + (1.0 - psi) * v2uc[node2] + psi * v2uc[node3]) / hz;
		}

		if(UseGradDivA)
		{
			v[0]+=(-(1.0-eta)*v2as[node0]+(1.0-eta)*v2as[node1]-eta*v2as[node2]+eta*v2as[node3])/hr;
			v[1]+=(-(1.0-eta)*v2ac[node0]+(1.0-eta)*v2ac[node1]-eta*v2ac[node2]+eta*v2ac[node3])/hr;
		}

		sig=sigma[rect[i].mtr-1];
		weps=w*dpr[rect[i].mtr-1];

		znam=sig*sig+weps*weps;

		rs=(sig*v[0]+weps*v[1])/znam;
		rc=(-weps*v[0]+sig*v[1])/znam;
	}
}

void SmoothOutput2D::GetSmoothVarForZR_U(int k, int jj, double *v2us, double *v2uc, double *v2as, double *v2ac, double &rs, double &rc)
{
	int i, node0, node1, node2, node3;
	double r, hr, hz, psi, eta;

	i = OutRecList[k].orobj.lem[jj];

	node0 = rect[i].nodes[0] - 1;
	node1 = rect[i].nodes[1] - 1;
	node2 = rect[i].nodes[2] - 1;
	node3 = rect[i].nodes[3] - 1;

	hr = pnt[node3].r - pnt[node0].r;
	hz = pnt[node3].z - pnt[node0].z;

	psi = 0.5;
	eta = (rec[k].z - pnt[node0].z) / hz;

	r = 0.5 * (pnt[node0].r + pnt[node3].r);

	rs = rc = 0.0;

	if (UseRotU)
	{
		rs += ((1.0 - psi) * (1.0 - eta) * v2us[node0] + psi * (1.0 - eta) * v2us[node1] + (1.0 - psi) * eta * v2us[node2] + psi * eta * v2us[node3]) / r +
			(-(1.0 - eta) * v2us[node0] + (1.0 - eta) * v2us[node1] - eta * v2us[node2] + eta * v2us[node3]) / hr;
		rc += ((1.0 - psi) * (1.0 - eta) * v2uc[node0] + psi * (1.0 - eta) * v2uc[node1] + (1.0 - psi) * eta * v2uc[node2] + psi * eta * v2uc[node3]) / r +
			(-(1.0 - eta) * v2uc[node0] + (1.0 - eta) * v2uc[node1] - eta * v2uc[node2] + eta * v2uc[node3]) / hr;
	}

	if (UseGradDivA)
	{
		rs += (-(1.0 - psi) * v2as[node0] - psi * v2as[node1] + (1.0 - psi) * v2as[node2] + psi * v2as[node3]) / hz;
		rc += (-(1.0 - psi) * v2ac[node0] - psi * v2ac[node1] + (1.0 - psi) * v2ac[node2] + psi * v2ac[node3]) / hz;
	}
}

double SmoothOutput2D::GetSmoothDzForZ(int k,int jj,double *v2,int OutType)
{
	int i,j,m,it,nn[2];
	double s,ls,v,hr,hz,psi;

	if(OutType && OutRecList[k].orfz[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].orfz[jj].ves) || (it && OutRecList[k].orfz[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].orfz[jj].lem[j+it];
					
					nn[0]=rect[m].nodes[0]-1;
					nn[1]=rect[m].nodes[3]-1;

					hr=pnt[nn[1]].r-pnt[nn[0]].r;
					hz=pnt[nn[1]].z-pnt[nn[0]].z;
					psi=0.5;


					v=(-(1.0-psi)*v2[rect[m].nodes[0]-1]-psi*v2[rect[m].nodes[1]-1]
						+(1.0-psi)*v2[rect[m].nodes[2]-1]+psi*v2[rect[m].nodes[3]-1])/hz;

					ls+=(OutRecList[k].orfz[jj].v[it][j]*v);
				}
				s+=ls*((it)? 1.0-OutRecList[k].orfz[jj].ves : OutRecList[k].orfz[jj].ves);
			}
		}
	}
	else
	{
		int node0,node1,node2,node3;

		i=OutRecList[k].oz.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;

		hr = pnt[node3].r - pnt[node0].r;
		hz = pnt[node3].z - pnt[node0].z;

		psi=(rec[k].r-pnt[node0].r)/hr;

		s=(-(1.0-psi)*v2[rect[i].nodes[0]-1]-psi*v2[rect[i].nodes[1]-1]
			+(1.0-psi)*v2[rect[i].nodes[2]-1]+psi*v2[rect[i].nodes[3]-1])/hz;
	}

	return s;
}

double SmoothOutput2D::GetSmoothDrForZ(int k,int jj,double *v2,int OutType)
{
	int i,j,m,it,nn[2];
	double s,ls,v,hr,hz,eta;

	if(OutType && OutRecList[k].orfz[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].orfz[jj].ves) || (it && OutRecList[k].orfz[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].orfz[jj].lem[j+it];
					
					nn[0]=rect[m].nodes[0]-1;
					nn[1]=rect[m].nodes[3]-1;

					hr = pnt[nn[1]].r - pnt[nn[0]].r;
					hz = pnt[nn[1]].z - pnt[nn[0]].z;
					
					eta=0.5;

					v=(-(1.0-eta)*v2[rect[m].nodes[0]-1]+(1.0-eta)*v2[rect[m].nodes[1]-1]
						-eta*v2[rect[m].nodes[2]-1]+eta*v2[rect[m].nodes[3]-1])/hr;
					ls+=(OutRecList[k].orfz[jj].v[it][j]*v);
				}
				s+=ls*((it)? 1.0-OutRecList[k].orfz[jj].ves : OutRecList[k].orfz[jj].ves);
			}
		}
	}
	else
	{
		int node0,node1,node2,node3;

		i=OutRecList[k].oz.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;

		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;

		eta=0.5;

		s=(-(1.0-eta)*v2[rect[i].nodes[0]-1]+(1.0-eta)*v2[rect[i].nodes[1]-1]
			-eta*v2[rect[i].nodes[2]-1]+eta*v2[rect[i].nodes[3]-1])/hr;
	}

	return s;
}

double SmoothOutput2D::GetSmoothValForB_U(int k,int jj,double *v2u,int OutType)
{
	int i,j,m,it,node0,node1,node2,node3;
	double s,ls,v,hr,hz,psi,eta,r;

	if(OutType || !OutRecList[k].orfz[jj].fl)
	{
		s=0.0;
		for(it=0;it<2;it++)
		{
			if((!it && OutRecList[k].orfz[jj].ves) || (it && OutRecList[k].orfz[jj].ves!=1.0))
			{
				ls=0.0;
				for(j=0;j<3;j++)
				{
					m=OutRecList[k].orfz[jj].lem[j+it];

					node0=rect[m].nodes[0]-1;
					node1=rect[m].nodes[1]-1;
					node2=rect[m].nodes[2]-1;
					node3=rect[m].nodes[3]-1;

					hr=pnt[node3].r-pnt[node0].r;
					hz=pnt[node3].z-pnt[node0].z;

					psi=eta=0.5;
					r=0.5*(pnt[node3].r+pnt[node0].r);

					v=0.0;
					
					v+=((1.0-psi)*(1.0-eta)*v2u[node0]+psi*(1.0-eta)*v2u[node1]+(1.0-psi)*eta*v2u[node2]+psi*eta*v2u[node3]);

					ls+=OutRecList[k].orfz[jj].v[it][j]*v;
				}
				s+=ls*((it)? 1.0-OutRecList[k].orfz[jj].ves : OutRecList[k].orfz[jj].ves);
			}
		}
	}
	else
	{
		s=0.0;
		i=OutRecList[k].oz.lem[jj];

		node0=rect[i].nodes[0]-1;
		node1=rect[i].nodes[1]-1;
		node2=rect[i].nodes[2]-1;
		node3=rect[i].nodes[3]-1;

		hr=pnt[node3].r-pnt[node0].r;
		hz=pnt[node3].z-pnt[node0].z;

		psi=(rec[k].r-pnt[node0].r)/hr;
		eta=0.5;
		r=rec[k].r;

		s+=((1.0-psi)*(1.0-eta)*v2u[node0]+psi*(1.0-eta)*v2u[node1]+(1.0-psi)*eta*v2u[node2]+psi*eta*v2u[node3]);
	}

	return s;
}
