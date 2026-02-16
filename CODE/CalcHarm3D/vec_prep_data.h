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
 *  This file contains headers for reading and storing edge mesh in 3D VFEM
 *  
 *  
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/vec_prep_data.h         
 *  Modified by D.Sc. Denis V. Vagin                                                                                 
 *  Novosibirsk State Technical University,                                                                          
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                 
 *  Version 2.0 December 10, 2024                                                                                    
*/

#pragma once
#include "T_Mapping.h"
#include "ElemNeib.h"

struct Point3D
{
	double x,y,z;
};

struct loc_source
{
	int elem;
	double x,y,z;
	double dx,dy,dz;
	double ksi,eta,dzeta;
	double d_ksi,d_eta,d_dzeta;
	double cur,len;
	loc_source(){elem=-1;}
};

struct LineSource
{
	Point3D A,B;
};

struct Tensor
{
	double val[3][3];

	Tensor()
	{
		Clear();
	}

	void Clear()
	{
		for(int i=0; i<3; i++)
		{
			for(int j=0; j<3; j++)
				val[i][j] = 0;
		}
	}

	Tensor& operator = (const double& d)
	{
		Clear();

		for(int j=0; j<3; j++)
			val[j][j] = d;

		return *this;
	}

	bool operator == (const double& d)
	{
		bool flag=true;
		int i, j;
		const double eps = 1e-6;

		for (i=0; i<3; i++)
		{
			for (j=0; j<3; j++)
			{
				if (i==j)
				{
					if(fabs(val[i][i] - d) > eps)
					{
						flag = false;
						break;
					}
				}
				else
				{
					if (val[i][j] != 0)
					{
						flag = false;
						break;
					}
				}
			}

			if (!flag)
				break;
		}

		return flag;
	}

	bool operator != (const double& d)
	{
		return !(*this == d);
	}
};

struct SigmaTensorRotate
{
	Tensor t;
	double m[3][3];
	SigmaTensorRotate()
	{
		int ii,jj;
		for (ii=0; ii<3; ii++)
		{
			for (jj=0; jj<3; jj++)
			{
				m[ii][jj]=0.0;
				t.val[ii][jj]=0.0;
			}
			m[ii][ii]=1.0;
		}
	}
	void GetMatrix(double r[3][3])
	{
		int ii,jj;
		for (ii=0; ii<3; ii++)
		{
			for (jj=0; jj<3; jj++)
			{
				r[ii][jj]=m[ii][jj];
			}
		}
	}
};

class Vec_Prep_Data
{
public:

	struct EnLine
	{
		long mtr;
		double es[3];
		double ec[3];
	};

	Vec_Prep_Data();
	~Vec_Prep_Data();

	int Read_prep_data();
	int ReadPrepDataHarmLoop(char *pointres_fname);
	int Read_mtz_1d();

	int Read_mesh_for_nonstat_problem(char *pointres_fname);
	int Read_infite0();

	int AnomalType,fda;

	int Read_3dmeshregular(long interval); 

	int LoadReceivers(char *fname, int& _n_pointres, double (*(&_pointres))[3]);

	long maxiter;
	double eps;
	long n_materials;
	double *mu3d;
	double *mu0;
	int n_pointresB;
	int n_pointresE;
	double (*pointresB)[3];
	double (*pointresE)[3];
	double *sigma3d;
	Tensor *sigma3dTensor,*sigma2dTensor;
	double *sigma0;
	SigmaTensorRotate *vSTR;
	double *dpr3d;
	double *dpr0;
	long kuzlov;
	long kpar;
	long kt1;
	long *l13d;
	long (*nver)[14];
	long *nvkat;
	double (*xyz)[3];
	long n_layers_1d; 
	double *layers_1d;
	double *sigma_1d;

	int LoadVectorE0ForLine(int n,int npls);
	vector< vector<EnLine> > EnForLine;

	int N_X,N_Y,N_Z;
	vector<int> regular;
	vector<double> Xcrd,Ycrd,Zcrd;

	long n_mesh_regular_x;
	long n_mesh_regular_y;
	double *mesh_regular_x;
	double *mesh_regular_y;

	long ntime;
	double *time;

 	double nu;
 	long alfa;
 	double *usin;
 	double *ucos;
 	double *z_1d;
 	long n_1d;

	T_Mapping_Vec *tmap;
	long tasktype;

	int npr, nfreq;

	double (*xyzt)[3];

	bool AnomalExists;
	vector<int> MaterialExists;

	int fdirect,nSrsRazb;
	vector<LineSource> lsrs;
	vector<vector<loc_source>> srs;

	vector<ElemNeib> ElemNeibVec;
};
