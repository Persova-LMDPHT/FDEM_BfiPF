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
 *  This file contains code for calculate 3D EM secondary field
 *
 *  Written by D.Sc. Denis V. Vagin 
 *  Novosibirsk State Technical University,
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024
 * 
*/

#include "stdafx.h"
#include "bound_cond_vec_harm.h"
#include "t_global_slae.h"
#include "T_Mapping.h"
#include "T_Portrait.h"
#include "give_out_vec_mt.h"
#include "pardiso.h"
#include "CheckInHex.h"
#include "in_out.h"
#include "ControlOMP.h"
#include "pcocr.h"

ofstream logfile;
ControlOMP omp;

const int size_i=sizeof(int);
const int size_d=sizeof(double);

int EFS[6][4]={{4,5,8,10},{0,2,8,9},{0,1,4,6},{6,7,9,11},{1,3,10,11},{2,3,5,7}};

bool IsFileExist(char *fname)
{
	bool flag;
	ifstream inf;

	flag=false;
	inf.open(fname);
	if(inf)
	{
		flag=true;
		inf.close();
	}
	inf.clear();

	return flag;
}

extern int Sn[6][4];

void GetGlobalCoordinates(Point3D *HexPnt,double *lc,double *gc)
{
	double ksi,eta,phi;

	ksi=lc[0];
	eta=lc[1];
	phi=lc[2];

	gc[0]=	HexPnt[0].x*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].x*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].x*(1-ksi)*(eta)*(1-phi)+HexPnt[3].x*(ksi)*(eta)*(1-phi)+
			HexPnt[4].x*(1-ksi)*(1-eta)*(phi)+HexPnt[5].x*(ksi)*(1-eta)*(phi)+
			HexPnt[6].x*(1-ksi)*(eta)*(phi)+HexPnt[7].x*(ksi)*(eta)*(phi);
	gc[1]=	HexPnt[0].y*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].y*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].y*(1-ksi)*(eta)*(1-phi)+HexPnt[3].y*(ksi)*(eta)*(1-phi)+
			HexPnt[4].y*(1-ksi)*(1-eta)*(phi)+HexPnt[5].y*(ksi)*(1-eta)*(phi)+
			HexPnt[6].y*(1-ksi)*(eta)*(phi)+HexPnt[7].y*(ksi)*(eta)*(phi);
	gc[2]=	HexPnt[0].z*(1-ksi)*(1-eta)*(1-phi)+HexPnt[1].z*(ksi)*(1-eta)*(1-phi)+
			HexPnt[2].z*(1-ksi)*(eta)*(1-phi)+HexPnt[3].z*(ksi)*(eta)*(1-phi)+
			HexPnt[4].z*(1-ksi)*(1-eta)*(phi)+HexPnt[5].z*(ksi)*(1-eta)*(phi)+
			HexPnt[6].z*(1-ksi)*(eta)*(phi)+HexPnt[7].z*(ksi)*(eta)*(phi);

}

void FindLocalCoordinates(Point3D &R,Point3D *HexPnt,double *lc)
{
	double EpsForFindLocalCoord = 1e-2;
	double CoeffForDiv = 1.2;
	int MaxDeep = 10;

	//logfile<<"Finding Reciver in Hehahedrob:"<<endl;

	int p,t,m,deep,crd[3][2]={{1,2},{0,2},{0,1}};
	double LocalCoord[2][3],CentGlob[3],CentLoc[3],DopLoc[3];
	double dist,disc,h[3];
	const double ods=0.16666666666666666;
	const double hds=0.33333333333333333;
	deep=0;
	LocalCoord[0][0]=0.0;LocalCoord[0][1]=0.0;LocalCoord[0][2]=0.0;
	LocalCoord[1][0]=1.0;LocalCoord[1][1]=1.0;LocalCoord[1][2]=1.0;
	CentLoc[0]=0.5;CentLoc[1]=0.5;CentLoc[2]=0.5;
	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
	disc=sqrt(
		(R.x-CentGlob[0])*(R.x-CentGlob[0])+
		(R.y-CentGlob[1])*(R.y-CentGlob[1])+
		(R.z-CentGlob[2])*(R.z-CentGlob[2])
		);
	h[0]=LocalCoord[1][0]-LocalCoord[0][0];
	h[1]=LocalCoord[1][1]-LocalCoord[0][1];
	h[2]=LocalCoord[1][2]-LocalCoord[0][2];
	do{
		if(disc<EpsForFindLocalCoord)break;

		for(m=0;m<3;m++){

			CentLoc[m]=LocalCoord[0][m]+ods*h[m];
			
			GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
			dist=sqrt(
					(R.x-CentGlob[0])*(R.x-CentGlob[0])+
					(R.y-CentGlob[1])*(R.y-CentGlob[1])+
					(R.z-CentGlob[2])*(R.z-CentGlob[2])
					);

			if(disc<dist){
				CentLoc[m]=LocalCoord[1][m]-ods*h[m];
			
				GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);
				dist=sqrt(
						(R.x-CentGlob[0])*(R.x-CentGlob[0])+
						(R.y-CentGlob[1])*(R.y-CentGlob[1])+
						(R.z-CentGlob[2])*(R.z-CentGlob[2])
						);
				if(dist<disc){
					disc=dist;
					LocalCoord[0][m]=LocalCoord[1][m]-hds*h[m];
				}
				else{
					LocalCoord[0][m]=LocalCoord[0][m]+hds*h[m];
					LocalCoord[1][m]=LocalCoord[1][m]-hds*h[m];
				}
			}
			else{
				disc=dist;
				LocalCoord[1][m]=LocalCoord[0][m]+hds*h[m];				
			}
			h[m]=LocalCoord[1][m]-LocalCoord[0][m];
			CentLoc[m]=0.5*(LocalCoord[0][m]+LocalCoord[1][m]);
		}
		deep++;
	}while(deep<MaxDeep);

	if(deep==MaxDeep){
		DopLoc[0]=CentLoc[0];
		DopLoc[1]=CentLoc[1];
		DopLoc[2]=CentLoc[2];
		do{
			t=0;
			for(m=0;m<3;m++){
				do{
					p=0;
					DopLoc[m]=CentLoc[m]-h[m];				
					if(DopLoc[m]>0){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
				do{
					p=0;
					DopLoc[m]=CentLoc[m]+h[m];				
					if(DopLoc[m]<1){
						GetGlobalCoordinates(HexPnt,DopLoc,CentGlob);
						dist=sqrt(
								(R.x-CentGlob[0])*(R.x-CentGlob[0])+
								(R.y-CentGlob[1])*(R.y-CentGlob[1])+
								(R.z-CentGlob[2])*(R.z-CentGlob[2])
								);
						if(dist<disc){
							disc=dist;
							CentLoc[m]=DopLoc[m];
							t=p=1;
						}
					}
				}while(p);
			}
		}while(t);
	}

	GetGlobalCoordinates(HexPnt,CentLoc,CentGlob);

	if(disc>EpsForFindLocalCoord){
		logfile<<"Reciver Eps Fail:"<<endl;
		logfile<<R.x<<'\t'<<R.y<<'\t'<<R.z<<endl;
		logfile<<CentGlob[0]<<'\t'<<CentGlob[1]<<'\t'<<CentGlob[2]<<endl;
		logfile<<sqrt(R.x*R.x+R.y*R.y+R.z*R.z)<<endl;
		logfile<<sqrt((R.x-CentGlob[0])*(R.x-CentGlob[0])+
				   (R.y-CentGlob[1])*(R.y-CentGlob[1])+
				   (R.z-CentGlob[2])*(R.z-CentGlob[2]))<<endl;
		logfile<<endl;
	}

	lc[0]=CentLoc[0];
	lc[1]=CentLoc[1];
	lc[2]=CentLoc[2];
}

int inputvec(long p_kpar,long &p_n,long &p_nc,long (**ed)[25],long (**edges)[2])
{
	int i,j;
	FILE *fp;
	ifstream inf;

	int p_vc,mv;

	inf.open("tsize3d_.dat");
	if(!inf){
		printf("Error open file tsize3d.dat");
		return 1;
	}
	inf>>p_vc;
	inf>>p_n;
	inf.close();
	inf.clear();

	p_nc=p_n-p_vc;

	if(!(fp=fopen("nodesforedges.dat","rb"))){
		printf("Error open file nodesforedges.dat");
		return 1;
	}

	if(!((*edges)=new long[p_n][2]))return 1;
	for(i=0;i<p_n;i++){
		fread((*edges)[i],size_i,2,fp);
		(*edges)[i][0]--;(*edges)[i][1]--;
	}
	fclose(fp);

	if(!(fp=fopen("edges.dat","rb"))){
		printf("Error open file edges.dat");
		return 1;
	}

	if(!((*ed)=new long[p_kpar][25]))return 1;
	for(i=0;i<p_kpar;i++){
		fread((*ed)[i],size_i,25,fp);
		for(j=0;j<25;j++)(*ed)[i][j]--;
	}
	fclose(fp);
	return 0;
}

void Unpuk(double *u,int npls,int slae_n,int tmap_n)
{
	int i,j,ipls;
	for(ipls=npls-1;ipls>0;ipls--)
	{
		for(i=0;i<slae_n;i++)
		{
			for(j=0;j<2;j++)
			{
				u[2*(ipls*tmap_n+slae_n-i-1)+j]=u[2*((ipls+1)*slae_n-i-1)+j];
			}
		}
	}
}

void normalize(Point3D &vec)
{
	double len=sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
	if(len>1e-6)
	{
		vec.x/=len;
		vec.y/=len;
		vec.z/=len;
	}
	else
	{
		vec.x=0.0;
		vec.y=0.0;
		vec.z=0.0;
	}
}

int GetNearestElement(loc_source &ls,Vec_Prep_Data *d)
{
	double scurr,sbest,m_s_best[8];
	int i,j,elem,ib,jb,i_best,n_best,m_i_best[8];
	Point3D Pmin,Pmax,Pcur,Hex[8],p,loc;
	if(!d->kpar){return -1;}
	p.x=ls.x;
	p.y=ls.y;
	p.z=ls.z;
	n_best=8;
	for(ib=0;ib<n_best;ib++)
	{
		m_s_best[ib]=1e+30;
		m_i_best[ib]=-1;
	}
	sbest=m_s_best[n_best-1];
	if(d->kpar<n_best){n_best=d->kpar;}

	for(elem=0;elem<d->kpar;elem++)
	{
		if(d->nvkat[elem])
		{
			j=d->nver[elem][0];
			Pcur.x=d->xyz[j][0];
			Pcur.y=d->xyz[j][1];
			Pcur.z=d->xyz[j][2];
			Hex[0]=Pmax=Pmin=Pcur;
			for(i=1;i<8;i++)
			{
				j=d->nver[elem][i];
				Pcur.x=d->xyz[j][0];
				Pcur.y=d->xyz[j][1];
				Pcur.z=d->xyz[j][2];
				Hex[i]=Pcur;
				if(Pcur.x<Pmin.x){Pmin.x=Pcur.x;}
				if(Pcur.y<Pmin.y){Pmin.y=Pcur.y;}
				if(Pcur.z<Pmin.z){Pmin.z=Pcur.z;}
				if(Pcur.x>Pmax.x){Pmax.x=Pcur.x;}
				if(Pcur.y>Pmax.y){Pmax.y=Pcur.y;}
				if(Pcur.z>Pmax.z){Pmax.z=Pcur.z;}
			}

			scurr=0.0;
			scurr += (p.x<Pmin.x)? (Pmin.x-p.x) : (p.x>Pmax.x)? (p.x-Pmax.x) : 0.0;
			scurr += (p.y<Pmin.y)? (Pmin.y-p.y) : (p.y>Pmax.y)? (p.y-Pmax.y) : 0.0;
			scurr += (p.z<Pmin.z)? (Pmin.z-p.z) : (p.z>Pmax.z)? (p.z-Pmax.z) : 0.0;

			if(scurr<sbest)
			{
				for(ib=0;ib<n_best;ib++)
				{
					if(m_i_best[ib]==-1 || scurr<m_s_best[ib])
					{
						break;
					}
				}

				if(m_i_best[ib]==-1)
				{
					m_i_best[ib]=elem;
					m_s_best[ib]=scurr;
				}
				else
				{
					for(jb=n_best-1;jb>ib;jb--)
					{
						m_i_best[jb]=m_i_best[jb-1];
						m_s_best[jb]=m_s_best[jb-1];
					}
					m_i_best[ib]=elem;
					m_s_best[ib]=scurr;
				}

				sbest=m_s_best[n_best-1];
			}
		}
	}

	i_best=-1;
	sbest=1e+30;
	for(ib=0;ib<n_best;ib++)
	{
		elem=m_i_best[ib];
		if(elem==-1){break;}
		for(i=0;i<8;i++)
		{
			j=d->nver[elem][i];
			Pcur.x=d->xyz[j][0];
			Pcur.y=d->xyz[j][1];
			Pcur.z=d->xyz[j][2];
			Hex[i]=Pcur;
		}
		loc.x = loc.y = loc.z = 0.5;
		if(!CheckInHex(Hex,p,loc))
		{
			scurr=0.0;
			scurr += (loc.x<0.0)? (-loc.x) : (loc.x>1.0)? (loc.x-1.0) : 0.0;
			scurr += (loc.y<0.0)? (-loc.y) : (loc.y>1.0)? (loc.y-1.0) : 0.0;
			scurr += (loc.z<0.0)? (-loc.z) : (loc.z>1.0)? (loc.z-1.0) : 0.0;
			if(scurr<sbest)
			{
				sbest=scurr;
				i_best=elem;
				ls.elem=elem;
				ls.ksi=loc.x;
				ls.eta=loc.y;
				ls.dzeta=loc.z;
			}
		}
	}

	return i_best;
}

int main()
{
	char str[256];
	int i,j,k,l,m;
	ifstream inf;
	int StartType,SolverType;
	ofstream ofp;
	double tmpd;
	pardiso_solver prds;

	double *h=NULL;
	double *y=NULL;
	double *y_omp=NULL;

	logfile.open("logharm3d");

	omp.InitNumThreads();

	StartType=0;
	inf.open("Harm3DParams");
	if(inf)
	{
		inf>>StartType;
		inf.close();
	}
	inf.clear();

	int p1,p2,ipls,npls,nprof;
	vector<int> RecvPlsIgB,RecvPlsIgE;

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
		inf>>j>>p1>>p2;
		npls+=p2-p1+1;
	}
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
		RecvPlsIgE[i+1]=RecvPlsIgE[i]+RecvPlsIgE[i+1];
		RecvPlsIgB[i+1]=RecvPlsIgB[i]+RecvPlsIgB[i+1];
	}

	Vec_Prep_Data *d=NULL;
	T_Mapping_Vec *tmap=NULL;
	
	if((d = new Vec_Prep_Data())==0) 
		throw logic_error("no memory for 3d mesh");
	
	d->tasktype=2;
	d->Read_prep_data();
	
	if((tmap = new T_Mapping_Vec(d->nver, d->xyz, d->kuzlov, d->kpar))==0)
	{
		if(d) {delete d; d=NULL;}
		throw logic_error("no memory for tmap");
	}

	if(inputvec(tmap->kpar,tmap->n,tmap->n_c,&(tmap->ed),&(tmap->edges)))
	{
		logfile<<"Error in inputvecmain"<<endl;
		logfile.close();
		return 1;
	}

	T_Portrait *p=NULL;
	VecBlockSLAE *slae=NULL;
	In_Out r;
	double *u=NULL;
	Bound_cond_vec_harm *bc=NULL;

	if(!(d->fdirect))
	{
		if (d->LoadVectorE0ForLine(tmap->n,npls)!=0)
		{
			cout<<"Can't read normal field"<<endl;
			return 1;
		}

		if(StartType==1)
		{
			ofp.open("e2.upk",ios::binary);
			for(ipls=0;ipls<npls;ipls++)
			{
				for(i=0;i<tmap->n;i++)
				{
					double es,ec;
					if((int)d->EnForLine[i+tmap->n*ipls].size())
					{
						es=d->EnForLine[i+tmap->n*ipls][0].es[0];
						ec=d->EnForLine[i+tmap->n*ipls][0].ec[0];
					}
					else
					{
						es=ec=0.0;
					}
					ofp.write((char *)&es,size_d);
					ofp.write((char *)&ec,size_d);
				}
			}
			ofp.close();
			ofp.clear();
		}

		if(!d->AnomalExists)
		{
			ofstream outf;
			double tmpd;

			tmpd=0.0;
			ofp.open("v3.dat",ios::binary);
			for(i=0;i<2*tmap->n_c*npls;i++){ofp.write((char *)&tmpd,size_d);}
			ofp.close();
			ofp.clear();
			ofp.open("v3.upk",ios::binary);
			for(i=0;i<2*tmap->n*npls;i++){ofp.write((char *)&tmpd,size_d);}
			ofp.close();
			ofp.clear();

			for(ipls=0;ipls<npls;ipls++)
			{
				sprintf(str,"e3d_anom.%d",ipls+1);
				outf.open(str);
				outf<<setiosflags(ios_base::scientific)<<setprecision(14);
				for(i=RecvPlsIgE[ipls];i<RecvPlsIgE[ipls+1];i++)
				{
					outf<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<'\n';
				}
				outf.close();
				outf.clear();

				sprintf(str,"b3d_anom.%d",ipls+1);
				outf.open(str);
				outf<<setiosflags(ios_base::scientific)<<setprecision(14);
				for(i=RecvPlsIgB[ipls];i<RecvPlsIgB[ipls+1];i++)
				{
					outf<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<' '<<0.0<<'\n';
				}
				outf.close();
				outf.clear();
			}

			return 0;
		}
	}
	else
	{
		if(StartType==1)
		{
			ofp.open("e2.upk",ios::binary);
			for(ipls=0;ipls<npls;ipls++)
			{
				for(i=0;i<tmap->n;i++)
				{
					double es,ec;
					es=ec=0.0;
					ofp.write((char *)&es,size_d);
					ofp.write((char *)&ec,size_d);
				}
			}
			ofp.close();
			ofp.clear();
		}

		vector<LineSource> &lsrs=d->lsrs;
		vector<vector<loc_source>> &srs=d->srs;
		int nSrsRazb=d->nSrsRazb;

		srs.resize(npls);

		inf.open("sours");
		if(!inf)
		{
			logfile<<"Error in open file "<<"sours"<<endl;
			cout<<"Error in open file "<<"sours"<<endl;
			return 1;
		}
		lsrs.resize(npls);
		for(i=0;i<npls;i++)
		{
			LineSource &ls=lsrs[i];
			inf>>ls.A.x>>ls.A.y>>ls.A.z;
			inf>>ls.B.x>>ls.B.y>>ls.B.z;
			srs[i].resize(nSrsRazb);
		}
		inf.close();
		inf.clear();

		int elem,nn;
		const double eps_loc_crd=1e-2;
		const double bnd_loc_min=0.0-eps_loc_crd;
		const double bnd_loc_max=1.0+eps_loc_crd;

		for(i=0;i<npls;i++)
		{
			LineSource &gen=lsrs[i];
			k=(int)srs[i].size();
			for(j=0;j<k;j++)
			{
				loc_source &ls=srs[i][j];
				Point3D dir,pp,pm;
				pm=gen.A;
				pp=gen.B;
				dir.x=pp.x-pm.x;
				dir.y=pp.y-pm.y;
				dir.z=pp.z-pm.z;
				dir.x/=nSrsRazb;
				dir.y/=nSrsRazb;
				dir.z/=nSrsRazb;
				ls.len=sqrt(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z);
				ls.x=pm.x+(j+0.5)*dir.x;
				ls.y=pm.y+(j+0.5)*dir.y;
				ls.z=pm.z+(j+0.5)*dir.z;
				ls.dx=dir.x;
				ls.dy=dir.y;
				ls.dz=dir.z;
			}
		}

		for(elem=0;elem<d->kpar;elem++)
		{
			if(d->nvkat[elem])
			{
				Point3D Pmin,Pmax,Pcur,Hex[8];

				j=d->nver[elem][0];
				Pcur.x=d->xyz[j][0];
				Pcur.y=d->xyz[j][1];
				Pcur.z=d->xyz[j][2];
				Hex[0]=Pmax=Pmin=Pcur;
				for(i=1;i<8;i++)
				{
					j=d->nver[elem][i];
					Pcur.x=d->xyz[j][0];
					Pcur.y=d->xyz[j][1];
					Pcur.z=d->xyz[j][2];
					Hex[i]=Pcur;
					if(Pcur.x<Pmin.x){Pmin.x=Pcur.x;}
					if(Pcur.y<Pmin.y){Pmin.y=Pcur.y;}
					if(Pcur.z<Pmin.z){Pmin.z=Pcur.z;}
					if(Pcur.x>Pmax.x){Pmax.x=Pcur.x;}
					if(Pcur.y>Pmax.y){Pmax.y=Pcur.y;}
					if(Pcur.z>Pmax.z){Pmax.z=Pcur.z;}
				}

				Pmin.x-=1e-3;
				Pmin.y-=1e-3;
				Pmin.z-=1e-3;
				Pmax.x+=1e-3;
				Pmax.y+=1e-3;
				Pmax.z+=1e-3;

				for(i=0;i<npls;i++)
				{
					LineSource &gen=lsrs[i];

					Point3D loc;

					loc.x=loc.y=loc.z=0.5;

					k=(int)srs[i].size();
					for(j=0;j<k;j++)
					{
						loc_source &ls=srs[i][j];
						if(ls.elem==-1)
						{
							const double delta=0.1;
							Point3D p,dir,pp,pm;

							pm=gen.A;
							pp=gen.B;
							p.x=ls.x;
							p.y=ls.y;
							p.z=ls.z;
							dir.x=ls.dx;
							dir.y=ls.dy;
							dir.z=ls.dz;

							if(p.x>=Pmin.x && p.x<=Pmax.x && p.y>=Pmin.y && p.y<=Pmax.y && p.z>=Pmin.z && p.z<=Pmax.z)
							{
								normalize(dir);
								pm.x=p.x-delta*dir.x;
								pm.y=p.y-delta*dir.y;
								pm.z=p.z-delta*dir.z;
								pp.x=p.x+delta*dir.x;
								pp.y=p.y+delta*dir.y;
								pp.z=p.z+delta*dir.z;

								nn=CheckInHex(Hex,p,loc);
								if (!nn &&
									(loc.x>bnd_loc_min && loc.x<bnd_loc_max &&
									loc.y>bnd_loc_min && loc.y<bnd_loc_max &&
									loc.z>bnd_loc_min && loc.z<bnd_loc_max ))
								{
									ls.elem=elem;
									ls.ksi=loc.x;
									ls.eta=loc.y;
									ls.dzeta=loc.z;

									nn=CheckInHex(Hex,pp,dir);
									if(nn)
									{
										logfile<<"Can't find local coordinates for point ";
										logfile<<pp.x<<' '<<pp.y<<' '<<pp.z<<' '<<endl;
										exit(1);
									}
									ls.d_ksi=dir.x;
									ls.d_eta=dir.y;
									ls.d_dzeta=dir.z;
									nn=CheckInHex(Hex,pm,dir);
									if(nn)
									{
										logfile<<"Can't find local coordinates for point ";
										logfile<<pm.x<<' '<<pm.y<<' '<<pm.z<<' '<<endl;
										exit(1);
									}
									ls.d_ksi-=dir.x;
									ls.d_eta-=dir.y;
									ls.d_dzeta-=dir.z;

									dir.x=ls.d_ksi;
									dir.y=ls.d_eta;
									dir.z=ls.d_dzeta;
									normalize(dir);
									ls.d_ksi=dir.x;
									ls.d_eta=dir.y;
									ls.d_dzeta=dir.z;
								}
							}
						}
					}
				}
			}
		}

		for(i=0;i<npls;i++)
		{
			k=(int)srs[i].size();
			for(j=0;j<k;j++)
			{
				loc_source &ls=srs[i][j];
				if(ls.elem==-1)
				{
					elem=GetNearestElement(ls,d);
					if(elem!=-1)
					{
						Point3D p,dir,pp,pm,Hex[8];
						const double delta=0.1;
						for(l=0;l<8;l++)
						{
							m=d->nver[elem][l];
							p.x=d->xyz[m][0];
							p.y=d->xyz[m][1];
							p.z=d->xyz[m][2];
							Hex[l]=p;
						}
						p.x=ls.x;
						p.y=ls.y;
						p.z=ls.z;
						dir.x=ls.dx;
						dir.y=ls.dy;
						dir.z=ls.dz;
						normalize(dir);
						pm.x=p.x-delta*dir.x;
						pm.y=p.y-delta*dir.y;
						pm.z=p.z-delta*dir.z;
						pp.x=p.x+delta*dir.x;
						pp.y=p.y+delta*dir.y;
						pp.z=p.z+delta*dir.z;

						nn=CheckInHex(Hex,pp,dir);
						if(nn)
						{
							logfile<<"Can't find local coordinates for point ";
							logfile<<pp.x<<' '<<pp.y<<' '<<pp.z<<' '<<endl;
							exit(1);
						}
						ls.d_ksi=dir.x;
						ls.d_eta=dir.y;
						ls.d_dzeta=dir.z;
						nn=CheckInHex(Hex,pm,dir);
						if(nn)
						{
							logfile<<"Can't find local coordinates for point ";
							logfile<<pm.x<<' '<<pm.y<<' '<<pm.z<<' '<<endl;
							exit(1);
						}
						ls.d_ksi-=dir.x;
						ls.d_eta-=dir.y;
						ls.d_dzeta-=dir.z;

						dir.x=ls.d_ksi;
						dir.y=ls.d_eta;
						dir.z=ls.d_dzeta;
						normalize(dir);
						ls.d_ksi=dir.x;
						ls.d_eta=dir.y;
						ls.d_dzeta=dir.z;
					}
					else
					{
						logfile<<"Error in find element for sours"<<endl;
						exit(1);
					}
				}
			}
		}
	}

	if((p = new T_Portrait((long*)tmap->ed, tmap->n, tmap->n_c, d->kpar))==0)
		Memory_allocation_error("p", "Loop_Harm_Vector_FEM");

	cout <<"p->Gen_Portrait();"<<endl;
	p->Gen_Portrait();

	cout <<"p->Gen_idi_ijg(d->nvkat);"<<endl;
	p->Gen_idi_ijg(d->nvkat, d->nver);

	if((slae = new VecBlockSLAE(p->ig, p->jg, p->idi, p->ijg, d->kpar, tmap->n, tmap->n_c,	
		d->xyz, d->nver, tmap->ed, tmap->edges, 
		d->nvkat, d->nu,  d->mu3d, d->mu0, d->sigma3d, d->sigma0, d->dpr3d, d->dpr0,
		d->alfa, d->n_1d, d->z_1d, d->usin, d->ucos, npls))==0)
		Memory_allocation_error("slae", "Loop_Harm_Vector_FEM");
	slae->tasktype=d->tasktype;
	cout <<"slae->AsmBlockSLAE();"<<endl;

	slae->AsmBlockSLAE_MATRIX(d);

	if((bc = new Bound_cond_vec_harm(d->kuzlov, tmap->n_c, d->kt1, d->l13d, tmap->edges, d->kpar))==0)
		Memory_allocation_error("bc", "Loop_Harm_Vector_FEM");

	bc->SetHomogenDirichletCond(p->ig, p->jg, p->idi, p->ijg, slae->di_block, slae->gg_block, slae->pr, npls);

	if((u = new double[2*(tmap->n)*npls])==0) Memory_allocation_error("u", "Loop_Harm_Vector_FEM");

	long n_last_iter;
	double eps_last_iter;
	double residual_decrease;
	double eps_x_change;
	bool error_last_iter;
	FILE *fp=NULL;

	double r_pr;

	error_last_iter = false;
	if((fp = fopen("stop_criteria","r"))==0)
	{
		error_last_iter = true;
	}
	else
	{
		int n;

		n = fscanf(fp, "%ld %lf %lf %lf", &n_last_iter, &eps_last_iter, &residual_decrease, &eps_x_change);

		if(n!=4)
			error_last_iter = true;

		if(n_last_iter<1)
			error_last_iter = true;

		fclose(fp);
	}
	if(error_last_iter)
	{
		eps_last_iter = 1.0001;
		n_last_iter = 30;
		residual_decrease = 1.1;
		eps_x_change = 1e-4;
	}

	SolverType = 0;
	inf.open("met_slae_3d");
	if (inf)
	{
		inf >> SolverType;
		inf >> d->eps;
		inf >> d->maxiter;
		inf.close();
	}
	inf.clear();

	if (!SolverType)
	{
		double eps_reg = 1e-12;

		inf.open("eps_reg");
		if (inf)
		{
			inf >> eps_reg;
			inf.close();
		}
		inf.clear();

		for (i = 0; i < slae->idi[slae->n_edges_c]; i++) { slae->di_block[i] += eps_reg * slae->di_block[i]; }
	}

	int nthrfrq,AnomalType,nthreads;
	
	AnomalType=0;
	nthrfrq=1;

	inf.open("AnomalType");
	if(inf)
	{
		inf>>AnomalType;
		inf.close();
	}
	inf.clear();

	if(!AnomalType)
	{
		inf.open("nthreads.txt");
		if(inf)
		{
			inf>>nthrfrq;
			inf.close();
		}
		inf.clear();
	}
	else
	{
		inf.open("nProcAnom.txt");
		if(inf)
		{
			inf>>nthrfrq;
			inf.close();
		}
		inf.clear();
	}

	omp.SetNumThreads(nthrfrq);

	slae->AsmBlockSLAE_PR(d);

	for(ipls=0;ipls<npls;ipls++)
	{
		bc->SetHomogenDirichletCondPr(p->ig, p->jg, p->idi, p->ijg, slae->di_block, slae->gg_block, slae->pr+2*bc->n_edges*ipls);
	}

	char fv3[256];
	sprintf(fv3, "v3.dat");

	inf.open(fv3,ios::binary);
	if(!inf)
	{
		inf.clear();

		double fnrm;

		fnrm=0.0;
		for(i=0;i<2*slae->n_edges_c*npls;i++){fnrm+=slae->pr[i]*slae->pr[i];}
		cout<<"fnrm= "<<fnrm<<endl;

		if(fnrm>1e-60)
		{
			if(!SolverType)
			{
				prds.factorize(slae->n_edges_c,(int *)p->ig,(int *)p->jg,slae->gg_block,slae->gg_block,slae->di_block, (int *)p->idi, (int *)p->ijg,nthrfrq);
				prds.solve_nrhs(npls,slae->pr,u);
			}
			else
			{
				double *sig = NULL;
				double *y_omp = NULL;
				long *ig_t = NULL;
				long *jg_t = NULL;
				double *gg_t = NULL;
				PCOCR pcocr;

				sig = new double[d->n_materials];
				y_omp = new double[slae->n_edges_c * 2 * omp.GetMaxThreads()];
				for (i = 0; i < d->n_materials; i++)
				{
					if (d->sigma3d[i] == 0)
						sig[i] = 0;
					else
						sig[i] = 1;
				}
				int n_nodes_c;
				int tsize3d;

				if(IsFileExist("tsize3d.dat"))
				{
					r.Read_Long_From_Txt_File("tsize3d.dat", (long*)&tsize3d);
				}
				else
				{
					tsize3d = 0;
				}
				n_nodes_c = d->kuzlov - tsize3d;

				ig_t = new long[d->kuzlov + 1];
				if (IsFileExist("ig3d.dat"))
				{
					r.Read_Bin_File_Of_Long("ig3d.dat", &ig_t[n_nodes_c], tsize3d + 1, 1);
					for (i = 0; i < n_nodes_c; i++) { ig_t[i] = 1; }
				}
				else
				{
					for (i = 0; i < d->kuzlov + 1; i++) { ig_t[i] = 1; }
				}
				for (i = 0; i < d->kuzlov + 1; i++){ig_t[i]--;}
				int ig_t_n_1 = ig_t[d->kuzlov];
				gg_t = new double[ig_t_n_1];
				if (IsFileExist("gg3d.dat"))
				{
					r.Read_Bin_File_Of_Double("gg3d.dat", gg_t, ig_t_n_1, 1);
				}
				jg_t = new long[ig_t_n_1];
				if (IsFileExist("jg3d.dat"))
				{
					r.Read_Bin_File_Of_Long("jg3d.dat", jg_t, ig_t_n_1, 1);
					for (i = 0; i < ig_t_n_1; i++) { jg_t[i]--; }
				}

				int* is_node_bound = NULL;

				if ((is_node_bound = new int[d->kuzlov]) == 0)
					Memory_allocation_error("is_node_bound", "Bound_cond_vec_harm::Make_list_of_bound_edges_harm");

				for (i = 0; i < d->kuzlov; i++)
				{
					is_node_bound[i] = 0;
				}

				for (i = 0; i < d->kt1; i++)
				{
					is_node_bound[d->l13d[i]] = 1;
				}

				for (i = 0; i < slae->n * npls; i++)
				{
					u[i] = 0.0;
				}

				for (ipls = 0; ipls < npls; ipls++)
				{
					double* _pr = slae->pr + ipls * slae->n;
					double* _u = u + ipls * slae->n;
					pcocr.PCOCR_2x2_Folded(slae->n_edges_c * 2, (int*)p->ig, (int*)p->jg, (int*)p->idi, (int*)p->ijg, slae->di_block, slae->gg_block, _pr, _u, d->eps, d->maxiter, y_omp, d->kpar, tmap->n_c, n_nodes_c, d->xyz, (int*)d->nvkat, (int(*)[14])d->nver, sig, (int(*)[2]) tmap->edges, (int*)ig_t, (int*)jg_t, gg_t, is_node_bound);
				}

				delete [] is_node_bound;
				delete [] sig;
				delete [] y_omp;
				delete [] ig_t;
				delete [] jg_t;
				delete [] gg_t;
			}
		}
		else
		{
			for(i=0;i<2*slae->n_edges_c*npls;i++){u[i]=0.0;}
		}

		r.Write_Bin_File_Of_Double((char*)(const char*)fv3, u, slae->n_edges_c*npls, 2);
	}
	else
	{
		inf.close();
		inf.clear();
		r.Read_Bin_File_Of_Double((char*)(const char*)fv3, u, slae->n_edges_c*npls, 2);
	}

	Unpuk(u,npls,tmap->n_c,tmap->n);
	
	r.Write_Bin_File_Of_Double("v3.upk", u, tmap->n*npls, 2);

	int nPntB,nPntE;

	nPntB=nPntE=0;

	Give_out_vec_mt *give_out_1=NULL;
	Give_out_vec_mt *give_out_2=NULL;

	Subdomain::npls=npls;

	double (*xyz0)[3];
	xyz0=NULL;
	d->xyzt=NULL;
	
	inf.open("xyz0.dat",ios::binary);
	if(inf)
	{
		xyz0=new double[d->kuzlov][3];
		for(i=0;i<d->kuzlov;i++)
		{
			inf>xyz0[i][0]>xyz0[i][1]>xyz0[i][2];
		}
		inf.close();
	}
	inf.clear();

	if(!xyz0)
	{
		inf.open("xyz.dat",ios::binary);
		if(inf)
		{
			xyz0=new double[d->kuzlov][3];
			for(i=0;i<d->kuzlov;i++)
			{
				inf>xyz0[i][0]>xyz0[i][1]>xyz0[i][2];
			}
			inf.close();
		}
	}

	d->LoadReceivers("xyzVectorB", d->n_pointresB, d->pointresB);
	d->n_pointresE=0;
	give_out_1 = new Give_out_vec_mt(d, tmap, u);
	if(d->n_pointresB)
	{
		nPntB=d->n_pointresB;

		d->xyzt=d->xyz;
		if(xyz0)
		{
			d->xyz=xyz0;
		}

		give_out_1->resultantB = new OutputResultant3d(give_out_1,vtWithoutDiscontinuity);
		give_out_1->resultantB->Prepare(6);

		if(xyz0)
		{
			d->xyz=d->xyzt;
			d->xyzt=NULL;
		}

		give_out_1->vvtb.resize(6);
		give_out_1->vvtb[0]=vtRotzASin;
		give_out_1->vvtb[1]=vtRotzACos;
		give_out_1->vvtb[2]=vtRotxASin;
		give_out_1->vvtb[3]=vtRotxACos;
		give_out_1->vvtb[4]=vtRotyASin;
		give_out_1->vvtb[5]=vtRotyACos;

		give_out_1->v3dat=u;
		for(ipls=0;ipls<npls;ipls++)
		{
			give_out_1->ipls_cur=ipls;
			give_out_1->Give_out_on_hex();
		}

		give_out_1->Write_B_to_files_for_harm_loop(StartType,RecvPlsIgB,0);
	}
	else
	{
		give_out_1->Write_NULL_B_to_files_for_harm_loop(0);
	}

	d->LoadReceivers("xyzVectorE", d->n_pointresE, d->pointresE);
	d->n_pointresB=0;
	give_out_2 = new Give_out_vec_mt(d, tmap, u);
	if(d->n_pointresE)
	{
		nPntE=d->n_pointresE;

		d->xyzt=d->xyz;
		if(xyz0)
		{
			d->xyz=xyz0;
		}

		give_out_2->resultantA = new OutputResultant3d(give_out_2,vtWithDiscontinuity);
		give_out_2->resultantA->Prepare(6);

		if(xyz0)
		{
			d->xyz=d->xyzt;
			d->xyzt=NULL;
		}

		give_out_2->vvta.resize(6);
		give_out_2->vvta[0]=vtAzSin;
		give_out_2->vvta[1]=vtAzCos;
		give_out_2->vvta[2]=vtAxSin;
		give_out_2->vvta[3]=vtAxCos;
		give_out_2->vvta[4]=vtAySin;
		give_out_2->vvta[5]=vtAyCos;

		give_out_2->v3dat=u;
		for(ipls=0;ipls<npls;ipls++)
		{
			give_out_2->ipls_cur=ipls;
			give_out_2->Give_out_on_hex();
		}

		give_out_2->Write_E_to_files_for_harm_loop(StartType,RecvPlsIgE,0);
	}
	else
	{
		give_out_2->Write_NULL_E_to_files_for_harm_loop(0);
	}

	if(nPntB)
	{
		give_out_1->resultantB->StopSolvers();
		delete give_out_1->resultantB;
		give_out_1->resultantB = NULL;
		if(give_out_1)
		{
			delete give_out_1;
			give_out_1=NULL;
		}
	}

	if(nPntE)
	{
		give_out_2->resultantA->StopSolvers();
		delete give_out_2->resultantA;
		give_out_2->resultantA = NULL;
		if(give_out_2)
		{
			delete give_out_2;
			give_out_2=NULL;
		}
	}

	if(xyz0)
	{
		delete [] xyz0;
		xyz0=NULL;
	}
	d->xyzt=NULL;

	if(!SolverType)
	{
		prds.stop_solver();
	}

	if(h) {delete [] h; h=NULL;}
	if(y) {delete [] y; y=NULL;}
	if(y_omp) {delete [] y_omp; y_omp=NULL;}

	if(p) {delete p; p=NULL;}
	if(u) {delete [] u; u=NULL;}
	if(bc) {delete bc; bc=NULL;}
	if(slae) {delete slae; slae=NULL;}
	if(d) {delete d; d=NULL;}
	if(tmap) {delete tmap; tmap=NULL;}

	logfile.close();
	logfile.clear();

	return 0;
}
