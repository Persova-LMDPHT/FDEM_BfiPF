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
 *  This file contains code of function for calculating the eleptic source part of field for harmonic task
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin 
 *  Novosibirsk State Technical University,                                             
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                   
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024                                                        
*/

#include "stdafx.h"
#include "Harm2D.h"
#include "rect_postproc_2d.h"
#include "Portret.h"
#include "in_out.h"
#include "global_slae_2d.h"
#include "For_Solvers.h"
#include "bound_cond_2d.h"
#include "iobinary.h"
#include "pardiso.h"

struct SqLoop
{
	double Ax,Ay,Az;
	double Bx,By,Bz;
};

int SolveHarmonic2D(Rect_preproc_data *d,double *v2,double cd,int npls, int nfreq)
{
	Portret *p=NULL;
	In_Out r;
	T_Global_SLAE_2d *a=NULL;
	vector<double> currentDeltaFunc;
	long k, j, i;
	double rd, zda, zdb;
	ifstream inf;

	int ipls,nthreads;
	vector<long> nodeDeltaFuncA,nodeDeltaFuncB;
	vector<SqLoop> GenSq;

	GenSq.resize(npls);
	nodeDeltaFuncA.resize(npls);
	nodeDeltaFuncB.resize(npls);
	currentDeltaFunc.resize(npls);

	nthreads=1;
	inf.open("nthreads.txt");
	if(inf)
	{
		inf>>nthreads;
		if(nthreads<1){nthreads=1;}
		inf.close();
	}
	inf.clear();

	inf.open("sours");
	if(!inf)
	{
		logfile<<"Error in open file "<<"sours"<<endl;
		cout<<"Error in open file "<<"sours"<<endl;
		return 1;
	}
	for(i=0;i<npls;i++){inf>>GenSq[i].Ax>>GenSq[i].Ay>>GenSq[i].Az>>GenSq[i].Bx>>GenSq[i].By>>GenSq[i].Bz;}
	inf.close();
	inf.clear();

	for(i=0;i<npls;i++){currentDeltaFunc[i]=1.0;}

	if((p = new Portret(d->nvtr, d->krect, d->kuzlov, d->n_nodes_c))==0)
		Memory_allocation_error("p","SolveHarmonicLoop2D");

	p->Gen_Portret_2d_rect_T();

	int *idi = NULL;
	idi = new int[d->n_nodes_c + 1];

	idi[0] = 0;
	for (i = 0; i < d->n_nodes_c; i++)
		idi[i + 1] = idi[i] + 2;

	int *ijg = NULL;
	ijg = new int[p->ig[d->n_nodes_c] + 1];
	ijg[0] = 0;
	for (i = 0; i < p->ig[d->n_nodes_c]; i++)
		ijg[i + 1] = ijg[i] + 2;

	for(ipls=0;ipls<npls;ipls++)
	{
		rd=0.0;
		zda=GenSq[ipls].Az;

		if (d->pType==Rect_preproc_data::ptNone||d->pType==Rect_preproc_data::ptAphi||d->pType==Rect_preproc_data::ptAphip)
				FindNodeWithDeltaFunction(d->kuzlov, d->xy, rd, zda, &(nodeDeltaFuncA[ipls]));
	}

	if (d->pType==Rect_preproc_data::ptHphi)
		for(i=0; i<d->n_materials; i++)
			swap(d->mu2d[i], d->sigma2d[i]);

	if((a = new T_Global_SLAE_2d(p->ig, p->jg, d->krect, d->kuzlov, d->kt1, d->nvtr, 
		d->xy, d->l1, d->nvk1, d->n_nodes_c, d->n_materials,
		d->mu2d, d->nvkat, d->nu, d->sigma2d, d->dpr2d, d->alfa, d->type_of_rect,
		d->sigma0, d->n_1d, d->coords_1d, d->sin_1d, d->cos_1d, cd, npls))==0)
		Memory_allocation_error("a","SolveHarmonicLoop2D");

	if(d->pType==Rect_preproc_data::ptHphi)a->type_of_task=1;else a->type_of_task=0;

	a->Assembling_With_T_Mapping(d->pType==Rect_preproc_data::ptAphi || d->pType==Rect_preproc_data::ptAphip);

	Bound_Cond *Bc=NULL;
	Bc = new Bound_Cond(a->di, a->ggl, a->ggu, a->pr, d->l1, d->nvk1, d->n_nodes_c, d->kt1, p->ig, p->jg);
	if(Bc == 0)
		Memory_allocation_error("Bc", "T_Global_SLAE_2d::Assembling_With_T_Mapping");
	Bc->Set_bound_cond_for_Hy_problem_with_symmetrization(npls, d->xy, cd, a->di_b, a->ggl_b, a->ggu_b);
	if(Bc) {delete Bc; Bc=NULL;}

	for(ipls=0;ipls<npls;ipls++)
	{
		if (d->pType==Rect_preproc_data::ptAphi||d->pType==Rect_preproc_data::ptAphip)
		{
			a->pr[(a->n*ipls+nodeDeltaFuncA[ipls])*2] = currentDeltaFunc[ipls]/(2.0*3.1415926535897932);
		}
	}

	char fv2[256];
	sprintf(fv2, "v2.dat_f%d", nfreq);

	for (i=0; i<d->n_nodes_c*npls; i++){v2[i*2]=v2[i*2+1]=0;}

	logfile << "nthreads= " << nthreads << endl;

	pardiso_solver prds;

	prds.factorize(d->n_nodes_c, (int*)p->ig, (int*)p->jg, a->ggl_b, a->ggu_b, a->di_b, idi, ijg, nthreads);
	prds.solve_nrhs(npls, a->pr, v2);
	prds.stop_solver();

	ofstream outf((const char*)fv2, ios::binary);
	ofstream outfs("v2s.dat", ios::binary);
	ofstream outfc("v2c.dat", ios::binary);

	for(ipls=0;ipls<npls;ipls++)
	{
		for (i=0; i<d->n_nodes_c; i++)
		{
				outf<v2[d->n_nodes_c*2*ipls+i*2];
				outfs<v2[d->n_nodes_c*2*ipls+i*2];
				outf<v2[d->n_nodes_c*2*ipls+i*2+1];
				outfc<v2[d->n_nodes_c*2*ipls+i*2+1];
		}
	}
	outf.close();
	outf.clear();
	outfs.close();
	outfs.clear();
	outfc.close();
	outfc.clear();

	if(p) {delete p; p=NULL;}
	if(a) {delete a; a=NULL;}

	return 0;	
}

 int ProcTask2DLine_Harmonic(double curr,int nfreq)
{
	int retc;

	Rect_preproc_data *d=NULL;
	double *v2=NULL;

	if ((d=new Rect_preproc_data)==0)
		throw logic_error("no memory");

	d->pType=Rect_preproc_data::ptAphi;

	d->Read_mesh_rz();
	
	ifstream inf;
	int i,k,p1,p2,npls,nprof;
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

	if ((v2=new double[d->n_nodes_c*2*npls])==0) 
		throw logic_error("no memory");

	if (retc=SolveHarmonic2D(d, v2, curr, npls, nfreq))
	{
		printf("Solving harmonic problem for line failed. Error code = %d.", retc);
		exit(1);
	}

	if (d) { delete d; d=NULL; }
	if (v2) { delete [] v2; v2=NULL; }

	return 0;
}

int FindNodeWithDeltaFunction(long kuzlov, double (*rz)[2], double rd, double zd, long *node)
{
	ifstream fin;
	double r, z;
	double r_cur, r_min;
	long i;

	*node = 0;
	r = rz[0][0];
	z = rz[0][1];
	r_cur = r_min = sqrt((r-rd)*(r-rd) + (z-zd)*(z-zd));
	for (i=1; i<kuzlov; i++)
	{
		r = rz[i][0];
		z = rz[i][1];
		r_cur = sqrt((r-rd)*(r-rd) + (z-zd)*(z-zd));

		if (r_cur < r_min)
		{
			r_min = r_cur;
			*node = i;
		}
	}

	return 0;
}
