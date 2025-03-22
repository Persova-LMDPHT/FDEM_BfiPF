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
 *  This file contains code of function for calculating the grounded source part of field for harmonic task
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                      
 *  Novosibirsk State Technical University,                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)     
 *  Version 2.0 December 10, 2024                                        
*/

#include "stdafx.h"
#include "MeshRZ.h"
#include "Portret.h"
#include "VelHarm2d.h"
#include "global_slae_2d.h"
#include "in_out.h"
#include "bound_cond_2d.h"
#include "For_Solvers.h"
#include "pardiso.h"

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

extern ofstream logfile;

int Output(MeshRZ &mesh, double *w_re, double *w_im,int npls)
{
	ifstream fin;
	ofstream ofp;
	int i, j;
	In_Out r;

	int nUall;
	double (*u_bnd_all)[2]=NULL;

	int *elmFlag=NULL;

	int nelmJ;
	int *elmJ=NULL;
	double (*sigJ)[2]=NULL;
	double (*epsJ)[2]=NULL;
	int *lowUpElm=NULL;

	int ipls,nthreads;

	nthreads=1;
	fin.open("nthreads.txt");
	if(fin)
	{
		fin>>nthreads;
		if(nthreads<1){nthreads=1;}
		fin.close();
	}
	fin.clear();

	double omega;
	r.Read_Double_From_Txt_File("nu", &omega);
	omega *= 2.0*PI;
	cout << "omega=" << omega << endl;

	fin.open("j_elm.txt");
	if (!fin)
	{
		Cannot_open_file("j_elm.txt", "");
	}
	fin >> nelmJ;
	fin.close();
	fin.clear();

	sigJ = new double[nelmJ][2];
	epsJ = new double[nelmJ][2];
	lowUpElm = new int[nelmJ];
	elmJ = new int[nelmJ];

	fin.open("jump.dat", ios_base::binary);
	if (!fin)
	{
		Cannot_open_file("jump.dat", "");
	}
	for(i=0;i<nelmJ;i++)
	{
		fin > elmJ[i];
		elmJ[i]--;
		fin > lowUpElm[i];
		lowUpElm[i]--;
		fin > sigJ[i][0] > epsJ[i][0];
		fin > sigJ[i][1] > epsJ[i][1];
		epsJ[i][0]*=omega;
		epsJ[i][1]*=omega;
	}
	fin.close();
	fin.clear();

	Portret p(mesh.nvtr, mesh.krect, mesh.kuzlov, mesh.n_nodes_c);
	p.Gen_Portret_2d_rect_T();

	int *idi = NULL;
	idi = new int[mesh.n_nodes_c+1];

	idi[0] = 0;
	for (i=0; i<mesh.n_nodes_c; i++)
		idi[i+1] = idi[i] + 2;

	int *ijg = NULL;
	ijg = new int[p.ig[mesh.n_nodes_c]+1];
	ijg[0] = 0;
	for (i=0; i<p.ig[mesh.n_nodes_c]; i++)
		ijg[i+1] = ijg[i] + 2;

	T_Global_SLAE_2d slae;
	slae.Asm_B_C_Output(p.ig, p.jg, mesh.krect, mesh.nvtr, mesh.xy,
		mesh.n_nodes_c, mesh.mu, mesh.nvkat, mesh.sigma, mesh.sigmaZ, mesh.dpr, omega, w_re, w_im, npls, mesh.kuzlov);

	int *is_node_bound=NULL;
	is_node_bound = new int[mesh.kuzlov];

	for (i=0; i<mesh.kuzlov; i++)
		is_node_bound[i] = 0;

	mesh.Find1kr(1e-4);
	for (i=0; i<mesh.kt1; i++)
		is_node_bound[mesh.l1[i]] = 1;

	nUall = 0;
	for (i=0; i<mesh.kuzlov; i++)
	{
		if (is_node_bound[i] != 0)
			nUall++;
	}

	u_bnd_all = new double[nUall][2];

	for (i=0; i<nUall; i++)
		u_bnd_all[i][0] = u_bnd_all[i][1] = 0;

	if (is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}

	double K=1e+4;
	Bound_Cond bc(slae.di_b,slae.ggl_b,slae.ggu_b,slae.pr,mesh.l1,mesh.n_nodes_c,mesh.kt1,p.ig,p.jg);
	
	bc.SetSigma(nelmJ,elmJ,sigJ,epsJ,lowUpElm,K,mesh.nvtr,mesh.xy,w_re,w_im,npls,mesh.kuzlov,mesh.nvkat,mesh.sigma,mesh.sigmaZ,mesh.dpr,omega);
	
	bc.Set_bound_cond_harm_output(u_bnd_all,npls);

	double *x=NULL;
	x = new double[mesh.n_nodes_c*2*npls];
	for(i=0;i<mesh.n_nodes_c*2*npls;i++){x[i]=0;}

	int n = mesh.n_nodes_c;

	pardiso_solver prds;

	cout<<"npls= "<<npls<<endl;

	prds.factorize(n, p.ig, p.jg, slae.ggl_b, slae.ggu_b, slae.di_b, idi, ijg, nthreads);
	prds.solve_nrhs(npls, slae.pr, x);
	prds.stop_solver();

	double *x_re=NULL;
	double *x_im=NULL;

	x_re=new double[mesh.kuzlov*npls];
	x_im=new double[mesh.kuzlov*npls];

	for(ipls=0;ipls<npls;ipls++)
	{
		for (int i=0; i<mesh.n_nodes_c; i++)
		{
			x_re[i+mesh.kuzlov*ipls] = x[i*2+mesh.n_nodes_c*2*ipls];
			x_im[i+mesh.kuzlov*ipls] = x[i*2+1+mesh.n_nodes_c*2*ipls];
		}
	}

	ofp.open("u.re",ios::binary);
	ofp.write((char *)x_re,mesh.kuzlov*npls*sizeof(double));
	ofp.close();
	ofp.clear();

	ofp.open("u.im",ios::binary);
	ofp.write((char *)x_im,mesh.kuzlov*npls*sizeof(double));
	ofp.close();
	ofp.clear();

	if (x_re) {delete [] x_re; x_re=NULL;}
	if (x_im) {delete [] x_im; x_im=NULL;}

	if (x) {delete [] x; x=NULL;}

	if (idi) {delete [] idi; idi=NULL;}
	if (ijg) {delete [] ijg; ijg=NULL;}

	if (elmJ) {delete [] elmJ; elmJ=NULL;}
	if (u_bnd_all) {delete [] u_bnd_all; u_bnd_all=NULL;}
	if (elmFlag) {delete [] elmFlag; elmFlag=NULL;}
	if (sigJ) {delete [] sigJ; sigJ=NULL;}
	if (lowUpElm) {delete [] lowUpElm; lowUpElm=NULL;}	

	return 0;
}

int OutputB()
{
	int ipls,npls,i,k,p1,p2,nprof;
	ifstream inf;

	logfile.open("log_output_b.txt");

	MeshRZ mesh;
	mesh.ReadMeshRZ();

	double *re = NULL;
	double *im = NULL; 

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

	re = new double[mesh.kuzlov*npls];
	im = new double[mesh.kuzlov*npls];

	inf.open("Ax\\v2s.dat",ios::binary);
	if(!inf)
	{
		logfile<<"Error in open file "<<"v2s.dat"<<endl;
		cout<<"Error in open file "<<"v2s.dat"<<endl;
		return 1;
	}
	inf.read((char *)re,mesh.kuzlov*npls*sizeof(double));
	inf.close();
	inf.clear();

	inf.open("Ax\\v2c.dat",ios::binary);
	if(!inf)
	{
		logfile<<"Error in open file "<<"v2c.dat"<<endl;
		cout<<"Error in open file "<<"v2c.dat"<<endl;
		return 1;
	}
	inf.read((char *)im,mesh.kuzlov*npls*sizeof(double));
	inf.close();
	inf.clear();

	Output(mesh,re,im,npls);

	if (re) {delete [] re; re=NULL;}
	if (im) {delete [] im; im=NULL;}

	return 0;
}
