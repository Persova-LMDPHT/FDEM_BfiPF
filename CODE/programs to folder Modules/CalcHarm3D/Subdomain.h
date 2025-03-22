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
 *  This file contains structures and headers of functions for building subdomains for smoothing 3D EM field
 *                                                                        
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/Subdomain.h
 *  Modified by D.Sc. Denis V. Vagin                                                                 
 *  Novosibirsk State Technical University,                                                          
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                 
 *  Version 2.0 December 10, 2024                                                                    
*/

#pragma once
#include "ElemNeib.h"
#include "Portret.h"
#include "T_Brick.h"
#include "AbstractFEM.h"

#include "pardiso.h"

struct _SIGMA_N{
	int N;
	int gtn,bf[4];
	double val[4];
};

struct Local_Vector
{
	double g[8];
	Local_Vector()
	{
		for(int i=0;i<8;i++)g[i]=0.0;
	}
};

class Subdomain
{
public:
	
	int my_ind,ncmp;
	char dir_name[32];

	static int subind;
	static int npls;
	static int ntimes;

	long material;

	vector<int> renumElemFromOldToNew;
	vector<int> renumElemFromNewToOld;
	vector<int> renumNodeFromNewToOld;
	vector< vector<int> > renumEdgeFromNewToOld;
	int n_elem;
	int n_nodes;
	int n_nodes_c;
	int n_nodes_dc;
	
	vector<double> ValueInCenter;	
	vector<ElemNeib> ElemNeibVec;

	long (*nver)[14];
	double (*xyz)[3];
	double (*xyz_r)[3];
	Portret *p;

	vector<int> ig_t,jg_t;
	vector<double> gg_t;

	void build_elem_neib_first(int nx,int ny,int nz,int *regular,vector<bool> &ElemCheck);
	void BuildSigmaStruct(vector<_SIGMA_N> &SIGMA_N,double (*xyz_p)[3]);
	void BuildTMatrix(vector<_SIGMA_N> &SIGMA_N);
	void calc_Tij(vector<_SIGMA_N> &SIGMA_N,int j,int li,vector<int> &iv,vector<double> &dv,int &t,double vklad);
	void CalcValuesAll_Node(double *v);
	void Correct_T_Node(int nc,vector<_SIGMA_N> &SIGMA_N,int i,int p,int *pnt,bool &Correct_T_Flag,double pVal1,double pVal2);

	double *di;
	double *gg;
	double *pr;
	double *x;

	double *d;
	double *sg;

	int Init(int material, int levelNeighbors, vector< vector<long> > &PointresForElem, AbstractFEM3D *TaskCalcMesh);

	void AsmGlobalMatrix();
	void AsmGlobalVector(int ipls,int icmp);
	void CalcRightPartVect(T_Brick &L);

	Subdomain();
	~Subdomain();

	vector<Local_Vector> vLV;

	pardiso_solver prds;

	int nthreads;
};
