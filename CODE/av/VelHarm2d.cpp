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
#include "pardiso.h"

extern ofstream logfile;

struct PointXYZ
{
	double x,y,z;
};

struct SqLoop
{
	PointXYZ A,B;
};

int read_harmonic_fields(int n,char *fs,char *fc,double *u)
{
	int i,m;
	ifstream inf;
	
	m=2*n;
	
	inf.open(fs,ios::binary);
	if(!inf)
	{
		cout<<"Error: can't open file "<<fs<<endl;
		logfile<<"Error: can't open file "<<fs<<endl;
		return 1;
	}
	for(i=0;i<m;i+=2){inf>u[i];}
	inf.close();
	inf.clear();

	inf.open(fc,ios::binary);
	if(!inf)
	{
		cout<<"Error: can't open file "<<fc<<endl;
		logfile<<"Error: can't open file "<<fc<<endl;
		return 1;
	}
	for(i=1;i<m;i+=2){inf>u[i];}
	inf.close();
	inf.clear();

	return 0;
}

int VelHarm2D_Div()
{
	int i,j;
	char f_name[20];
	In_Out r;
	int n,m,l;
	double eps;
	int maxiter;

	logfile.open("log_vel2d_harm.txt");

	int p1,p2,ipls,npls,nprof;
	ifstream inf;
	vector<SqLoop> GenSq;

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

	system("copy sigmaZ Ax\\sigmaZ");

	MeshRZ mesh;
	mesh.ReadMeshRZ("Ax");

	n = mesh.n_nodes_c*3;

	double *ar0harm = NULL;
	m=mesh.kuzlov*2*npls;
	ar0harm = new double[m];
	for (i=0; i<m; i++){ar0harm[i] = 0;}

	read_harmonic_fields(mesh.kuzlov*npls,"Ax\\v2s.dat","Ax\\v2c.dat",ar0harm);

	double omega;
	r.Read_Double_From_Txt_File("nu", &omega);
	omega *= 2.0*PI;
	cout << "omega=" << omega << endl;

	vector<double> curr;

	curr.resize(npls);

	inf.open("Ax\\currentval");
	if (!inf)
	{
		logfile << "Error in open file " << "Ax\\currentval" << endl;
		cout << "Error in open file " << "Ax\\currentval" << endl;
		return 1;
	}
	for (ipls = 0; ipls < npls; ipls++)
	{
		inf >> curr[ipls];
	}
	inf.close();
	inf.clear();

	Portret p(mesh.nvtr, mesh.krect, mesh.kuzlov, mesh.n_nodes_c);
	p.Gen_Portret_2d_rect_T_3x3();

	int *idi = NULL;
	idi = new int[n+1];

	idi[0] = 0;
	for (i=0; i<n; i++)
		idi[i+1] = idi[i] + 2;

	int *ijg = NULL;
	ijg = new int[p.ig[n]+1];
	ijg[0] = 0;
	for (i=0; i<p.ig[n]; i++)
		ijg[i+1] = ijg[i] + 2;

	T_Global_SLAE_2d slae;
	slae.AsmHarmDivMatrix(p.ig, p.jg, mesh.krect, mesh.nvtr, mesh.xy, mesh.n_nodes_c, mesh.mu, mesh.nvkat, mesh.sigma_xy, mesh.sigma_z, omega, ar0harm);
	slae.AsmHarmDivPr(p.ig, p.jg, mesh.krect, mesh.nvtr, mesh.xy, mesh.n_nodes_c, mesh.kuzlov,
		mesh.mu, mesh.nvkat, mesh.sigma_xy, mesh.sigma_z, omega, ar0harm, npls);

	Bound_Cond bc(slae.di_b, slae.ggl_b, slae.ggu_b, slae.pr, mesh.l1, mesh.n_nodes_c, mesh.kt1, p.ig, p.jg);

	for(ipls=0;ipls<npls;ipls++)
	{
		int ib;
		double bv,cv;
		ib=-1;
		bv=1e+30;
		for(i=0;i<mesh.qz;i++)
		{
			cv=fabs(GenSq[ipls].A.z-mesh.zm[i]);
			if(cv<bv)
			{
				ib=i;
				bv=cv;
			}
		}
		if(ib!=-1)
		{
			slae.pr[ib*mesh.qr*6+ipls*mesh.n_nodes_c*6+4] += curr[ipls]/(2.0*PI)/omega;
		}
	}

	bc.Set_bound_cond_harm(npls);
	bc.SwitchOffEquationsHarm(0,0,0);

	double *x=NULL;
	m=n*2*npls;
	x = new double[m];
	for (i=0; i<m; i++){x[i] = 0;}

	pardiso_solver prds;

	prds.factorize(n, p.ig, p.jg, slae.ggl_b, slae.ggu_b, slae.di_b, idi, ijg, 1);
	prds.solve_nrhs(npls, slae.pr, x);
	prds.stop_solver();

	double *ar_re=NULL;
	double *ar_im=NULL;
	double *az_re=NULL;
	double *az_im=NULL;
	double *v_re=NULL;
	double *v_im=NULL;

	m=mesh.kuzlov*npls;

	ar_re = new double[m];
	ar_im = new double[m];
	az_re = new double[m];
	az_im = new double[m];
	v_re = new double[m];
	v_im = new double[m];

	for(ipls=0;ipls<npls;ipls++)
	{
		m=n*2*ipls;
		l=mesh.kuzlov*ipls;

		for (i=0; i<mesh.n_nodes_c; i++)
		{
			ar_re[l+i] = x[m+i*6];
			ar_im[l+i] = x[m+i*6+1];
			az_re[l+i] = x[m+i*6+2];
			az_im[l+i] = x[m+i*6+3];
			v_re[l+i]  = x[m+i*6+4];
			v_im[l+i]  = x[m+i*6+5];
		}
	}

	m=mesh.kuzlov*npls;

	r.Write_Bin_File_Of_Double("ar.re", ar_re, m, 1);
	r.Write_Bin_File_Of_Double("ar.im", ar_im, m, 1);
	r.Write_Bin_File_Of_Double("az.re", az_re, m, 1);
	r.Write_Bin_File_Of_Double("az.im", az_im, m, 1);
	r.Write_Bin_File_Of_Double("v.re",   v_re, m, 1);
	r.Write_Bin_File_Of_Double("v.im",   v_im, m, 1);

	if (ar_re) {delete [] ar_re; ar_re=NULL;}
	if (ar_im) {delete [] ar_im; ar_im=NULL;}
	if (az_re) {delete [] az_re; az_re=NULL;}
	if (az_im) {delete [] az_im; az_im=NULL;}
	if (v_re) {delete [] v_re; v_re=NULL;}
	if (v_im) {delete [] v_im; v_im=NULL;}

	if (ar0harm) {delete [] ar0harm; ar0harm=NULL;}
	if (idi) {delete [] idi; idi=NULL;}
	if (ijg) {delete [] ijg; ijg=NULL;}
	if (x) {delete [] x; x=NULL;}

	logfile.close();
	logfile.clear();

	return 0;
}
