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
 *  This file contains code of functions calculating 2D local matirxes
 * 
 * 
 *  Written by  Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 December 10, 2024                                            
*/

#include "stdafx.h"
#include "Local_matrix_rz.h"

LocalMatrixRZ::LocalMatrixRZ(double r0, double r1, double z0, double z1, double mu, double sigma, double beta)
{
	this->r0 = rk = r0;
	this->r1 = r1;
	this->z0 = zk = z0;
	this->z1 = z1;

	hr = r1 - r0;
	hz = z1 - z0;

	if (hr <= 0 || hz <= 0)
	{
		cout << "hr<=0 or hz<=0 \n";
		exit(1);
	}

	this->mu = mu;
	this->sigma = sigma;
	this->beta = beta;
}
//------------------------------------------------------ 
// Calculate the local stiffness matrix of a finite element with a term from the rotor
//------------------------------------------------------ 
void LocalMatrixRZ::CalcLocalMatrixB_r2()
{
	GetStiffnessMatrix(r0, z1, r1, z0);
}
//------------------------------------------------------ 
// Calculate the local stiffness matrix of an anisotropic finite element
//------------------------------------------------------ 
void LocalMatrixRZ::CalcLocalMatrixB(double sig_xy, double sig_z)
{
	GetStiffnessMatrix0(r0, z1, r1, z0, sig_xy, sig_z);
	return;
}
//------------------------------------------------------ 
// Calculate the local stiffness matrix of a finite element  
//------------------------------------------------------       
void LocalMatrixRZ::CalcLocalMatrixB()
{
	GetStiffnessMatrix0(r0, z1, r1, z0);
	return;

	double t3 = 1.0/hr/hz;
	double t4 = hr*hr;
	double t5 = t4*hr;
	double t6 = rk*t4;
	double t7 = 4.0*t6;
	double t8 = hz*hz;
	double t9 = rk*t8;
	double t10 = 4.0*t9;
	double t11 = hr*t8;
	double t12 = 2.0*t11;
	double t15 = t3*(t5+t7+t10+t12)/12.0;
	double t16 = 2.0*t6;
	double t19 = t3*(t5+t16-t10-t12)/12.0;
	double t20 = 2.0*t9;
	double t23 = t3*(t5+t7-t20-t11)/12.0;
	double t26 = t3*(t5+t16+t20+t11)/12.0;
	double t27 = 3.0*t5;
	double t30 = t3*(t27+t7+t10+t12)/12.0;
	double t33 = t3*(t27+t7-t20-t11)/12.0;

	b[0][0] = t15;
	b[0][1] = t19;
	b[0][2] = -t23;
	b[0][3] = -t26;
	b[1][0] = t19;
	b[1][1] = t30;
	b[1][2] = -t26;
	b[1][3] = -t33;
	b[2][0] = -t23;
	b[2][1] = -t26;
	b[2][2] = t15;
	b[2][3] = t19;
	b[3][0] = -t26;
	b[3][1] = -t33;
	b[3][2] = t19;
	b[3][3] = t30;
}
//------------------------------------------------------ 
// Calculate the local mass matrix of the finite element
//------------------------------------------------------ 
void LocalMatrixRZ::CalcLocalMatrixC()
{
	GetMassMatrix(r0, z1, r1, z0);
	return;

	double t1 = hr*hz;
	double t2 = 4.0*rk;
	double t4 = t1*(t2+hr);
	double t5 = t4/36.0;
	double t8 = t1*(2.0*rk+hr);
	double t9 = t8/36.0;
	double t10 = t4/72.0;
	double t11 = t8/72.0;
	double t14 = t1*(t2+3.0*hr);
	double t15 = t14/36.0;
	double t16 = t14/72.0;

	c[0][0] = t5;
	c[0][1] = t9;
	c[0][2] = t10;
	c[0][3] = t11;
	c[1][0] = t9;
	c[1][1] = t15;
	c[1][2] = t11;
	c[1][3] = t16;
	c[2][0] = t10;
	c[2][1] = t11;
	c[2][2] = t5;
	c[2][3] = t9;
	c[3][0] = t11;
	c[3][1] = t16;
	c[3][2] = t9;
	c[3][3] = t15;
}
//------------------------------------------------------ 
// Calculate the local matrix from integrals with Psi*d_Psi/d_r of the finite element
//------------------------------------------------------ 
void LocalMatrixRZ::CalcLocalMatrixDr()
{
	double t1 = 3.0*rk;
	double t3 = hz*(t1+hr);
	double t4 = t3/18.0;
	double t5 = t3/36.0;
	double t8 = hz*(t1+2.0*hr);
	double t9 = t8/18.0;
	double t10 = t8/36.0;

	dr[0][0] = -t4;
	dr[0][1] = t4;
	dr[0][2] = -t5;
	dr[0][3] = t5;
	dr[1][0] = -t9;
	dr[1][1] = t9;
	dr[1][2] = -t10;
	dr[1][3] = t10;
	dr[2][0] = -t5;
	dr[2][1] = t5;
	dr[2][2] = -t4;
	dr[2][3] = t4;
	dr[3][0] = -t10;
	dr[3][1] = t10;
	dr[3][2] = -t9;
	dr[3][3] = t9;
}
//------------------------------------------------------ 
// Calculate the local matrix from the integrals with Psi*d_Psi/d_z of the finite element
//------------------------------------------------------ 
void LocalMatrixRZ::CalcLocalMatrixDz()
{
	double t1 = 4.0*rk;
	double t4 = hr*(t1+hr)/24.0;
	double t8 = hr*(2.0*rk+hr)/24.0;
	double t12 = hr*(t1+3.0*hr)/24.0;

	dz[0][0] = -t4;
	dz[0][1] = -t8;
	dz[0][2] = t4;
	dz[0][3] = t8;
	dz[1][0] = -t8;
	dz[1][1] = -t12;
	dz[1][2] = t8;
	dz[1][3] = t12;
	dz[2][0] = -t4;
	dz[2][1] = -t8;
	dz[2][2] = t4;
	dz[2][3] = t8;
	dz[3][0] = -t8;
	dz[3][1] = -t12;
	dz[3][2] = t8;
	dz[3][3] = t12;
}

void LocalMatrixRZ::CalcLocalMatrixHarm(double omega)
{
	CalcLocalMatrixB();
	CalcLocalMatrixB_r2();
	CalcLocalMatrixC();
	CalcLocalMatrixDr();
	CalcLocalMatrixDz();

	int i, j;

	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			ah[i*3  ][j*3  ][0] =  b_r2[i][j]/mu;
			ah[i*3  ][j*3+1][0] =  0.0;
			ah[i*3  ][j*3+2][0] =  dr[i][j]*sigma;

			ah[i*3+1][j*3  ][0] =  0.0;
			ah[i*3+1][j*3+1][0] =  b[i][j]/mu;
			ah[i*3+1][j*3+2][0] =  dz[i][j]*sigma;

			ah[i*3+2][j*3  ][0] =  0.0;
			ah[i*3+2][j*3+1][0] =  0.0;
			ah[i*3+2][j*3+2][0] =  b[i][j]*sigma/omega;


			ah[i*3  ][j*3  ][1] =  c[i][j]*sigma*omega;
			ah[i*3  ][j*3+1][1] =  0.0;
			ah[i*3  ][j*3+2][1] =  0.0;

			ah[i*3+1][j*3  ][1] =  0.0;
			ah[i*3+1][j*3+1][1] =  c[i][j]*sigma*omega;
			ah[i*3+1][j*3+2][1] =  0.0;

			ah[i*3+2][j*3  ][1] =  dr[j][i]*sigma;
			ah[i*3+2][j*3+1][1] =  dz[j][i]*sigma;
			ah[i*3+2][j*3+2][1] =  0.0;
		}
	}
}

void LocalMatrixRZ::CalcLocalMatrixHarm(double omega, double sig_xy, double sig_z)
{
	CalcLocalMatrixB();
	CalcLocalMatrixB_r2();
	CalcLocalMatrixB(sig_xy, sig_z);
	CalcLocalMatrixC();
	CalcLocalMatrixDr();
	CalcLocalMatrixDz();

	int i, j;

	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			ah[i*3  ][j*3  ][0] =  b_r2[i][j]/mu;
			ah[i*3  ][j*3+1][0] =  0.0;
			ah[i*3  ][j*3+2][0] =  dr[i][j]*sig_xy;

			ah[i*3+1][j*3  ][0] =  0.0;
			ah[i*3+1][j*3+1][0] =  b[i][j]/mu;
			ah[i*3+1][j*3+2][0] =  dz[i][j]*sig_z;

			ah[i*3+2][j*3  ][0] =  0.0;
			ah[i*3+2][j*3+1][0] =  0.0;
			ah[i*3+2][j*3+2][0] =  ba[i][j]/omega;

			ah[i*3  ][j*3  ][1] =  c[i][j]*sig_xy*omega;
			ah[i*3  ][j*3+1][1] =  0.0;
			ah[i*3  ][j*3+2][1] =  0.0;

			ah[i*3+1][j*3  ][1] =  0.0;
			ah[i*3+1][j*3+1][1] =  c[i][j]*sig_z*omega;
			ah[i*3+1][j*3+2][1] =  0.0;

			ah[i*3+2][j*3  ][1] =  dr[j][i]*sig_xy;
			ah[i*3+2][j*3+1][1] =  dz[j][i]*sig_z;
			ah[i*3+2][j*3+2][1] =  0.0;
		}
	}
}

void LocalMatrixRZ::CalcLocalMatrix()
{
	CalcLocalMatrixB();
	CalcLocalMatrixB_r2();
	CalcLocalMatrixC();
	CalcLocalMatrixDr();
	CalcLocalMatrixDz();

	int i, j;

	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			ab[i*3  ][j*3  ] =  b_r2[i][j]/mu;
			ab[i*3  ][j*3+1] =  0.0;
			ab[i*3  ][j*3+2] = dr[i][j]*sigma;

			ab[i*3+1][j*3  ] =  0.0;
			ab[i*3+1][j*3+1] =  b[i][j]/mu;
			ab[i*3+1][j*3+2] = dz[i][j]*sigma;

			ab[i*3+2][j*3  ] = 0.0;
			ab[i*3+2][j*3+1] = 0.0;
			ab[i*3+2][j*3+2] =  b[i][j]*sigma/beta;
		}
	}

	for (i=0; i<4; i++)
	{
		for (j=0; j<4; j++)
		{
			ac[i*3  ][j*3  ] =  c[i][j]*sigma;
			ac[i*3  ][j*3+1] =  0.0;
			ac[i*3  ][j*3+2] =  0.0;

			ac[i*3+1][j*3  ] =  0.0;
			ac[i*3+1][j*3+1] =  c[i][j]*sigma;
			ac[i*3+1][j*3+2] =  0.0;

			ac[i*3+2][j*3  ] = dr[j][i]*sigma/beta;
			ac[i*3+2][j*3+1] = dz[j][i]*sigma/beta;
			ac[i*3+2][j*3+2] =  0.0;
		}
	}

	return;
}

void LocalMatrixRZ::CalcLocalMatrix2(double beta)
{
	CalcLocalMatrix();
	int i, j;

	for (i=0; i<12; i++)
	{
		for(j=0; j<12; j++)
		{
			aa[i][j] = ab[i][j] + ac[i][j]*beta; 
		}
	}

	return;
}

void LocalMatrixRZ::CalcLocalVector(double *val_r, double *val_z, double *val_v, double beta)
{
	int i, j;
	double fr[4];
	double fz[4];
	double fv[4];
	double drt[4][4];

	for (i=0; i<4; i++)
		for(j=0; j<4; j++)
			drt[i][j] = dr[j][i];

	Mult_Plot((double*)c,   val_r, fr, 4);
	Mult_Plot((double*)c,   val_z, fz, 4);
	Mult_Plot((double*)drt, val_v, fv, 4);

	g[0] = fr[0];
	g[1] = fz[0];
	g[2] = fv[0]*sigma/beta;

	g[3] = fr[1];
	g[4] = fz[1];
	g[5] = fv[1]*sigma/beta;

	g[6] = fr[2];
	g[7] = fz[2];
	g[8] = fv[2]*sigma/beta;

	g[9]  = fr[3];
	g[10] = fz[3];
	g[11] = fv[3]*sigma/beta;
}

void LocalMatrixRZ::CalcLocalVectorM(double *val_r, double *val_z, double *val_v, double beta)
{
	int i, j;
	double fr[4];
	double fz[4];
	double fvr[4];
	double fvz[4];
	double drt[4][4];
	double dzt[4][4];

	for (i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			drt[i][j] = dr[j][i];
			dzt[i][j] = dz[j][i];
		}
	}

	Mult_Plot((double*)c,   val_r, fr, 4);
	Mult_Plot((double*)c,   val_z, fz, 4);
	Mult_Plot((double*)drt, val_r, fvr, 4);
	Mult_Plot((double*)dzt, val_z, fvz, 4);

	g[0] = fr[0]*sigma;
	g[1] = fz[0]*sigma;
	g[2] = (fvr[0] + fvz[0])*sigma/beta;

	g[3] = fr[1]*sigma;
	g[4] = fz[1]*sigma;
	g[5] = (fvr[1] + fvz[1])*sigma/beta;

	g[6] = fr[2]*sigma;
	g[7] = fz[2]*sigma;
	g[8] = (fvr[2] + fvz[2])*sigma/beta;

	g[9]  = fr[3]*sigma;
	g[10] = fz[3]*sigma;
	g[11] = (fvr[3] + fvz[3])*sigma/beta;
}

void LocalMatrixRZ::CalcLocalVectorHarm(double *val_r_re, double *val_z_re, double *val_v_re,
						 double *val_r_im, double *val_z_im, double *val_v_im)
{
	int i, j;

	double fr_re[4];
	double fz_re[4];
	double fvr_re[4];
	double fvz_re[4];

	double fr_im[4];
	double fz_im[4];
	double fvr_im[4];
	double fvz_im[4];

	double drt[4][4];
	double dzt[4][4];

	for (i=0; i<4; i++)
	{
		for(j=0; j<4; j++)
		{
			drt[i][j] = dr[j][i];
			dzt[i][j] = dz[j][i];
		}
	}

	Mult_Plot((double*)drt, val_r_re, fvr_re, 4);
	Mult_Plot((double*)drt, val_r_im, fvr_im, 4);

	gh[0][0] = 0;
	gh[1][0] = 0;
	gh[2][0] = (fvr_im[0])*sigma;

	gh[3][0] = 0;
	gh[4][0] = 0;
	gh[5][0] = (fvr_im[1])*sigma;

	gh[6][0] = 0;
	gh[7][0] = 0;
	gh[8][0] = (fvr_im[2])*sigma;

	gh[9][0]  = 0;
	gh[10][0] = 0;
	gh[11][0] = (fvr_im[3])*sigma;

	gh[0][1] = 0;
	gh[1][1] = 0;
	gh[2][1] = -(fvr_re[0])*sigma;

	gh[3][1] = 0;
	gh[4][1] = 0;
	gh[5][1] = -(fvr_re[1])*sigma;

	gh[6][1] = 0;
	gh[7][1] = 0;
	gh[8][1] = -(fvr_re[2])*sigma;

	gh[9][1]  = 0;
	gh[10][1] = 0;
	gh[11][1] = -(fvr_re[3])*sigma;
}

void LocalMatrixRZ::CalcLocalVectorHarmDiv(double *val_r_re, double *val_r_im, double sig_xy, double sig_z)
{
	double fvr_re[4];
	double fvr_im[4];

	Mult_Plot((double*)c, val_r_re, fvr_re, 4);
	Mult_Plot((double*)c, val_r_im, fvr_im, 4);

	gh[0][0] = 0; 
	gh[1][0] = 0; 
	gh[2][0] = fvr_im[0]*sig_xy;

	gh[3][0] = 0;
	gh[4][0] = 0;
	gh[5][0] = fvr_im[1]*sig_xy;

	gh[6][0] = 0;
	gh[7][0] = 0;
	gh[8][0] = fvr_im[2]*sig_xy;

	gh[9][0]  = 0;
	gh[10][0] = 0;
	gh[11][0] = fvr_im[3]*sig_xy;

	gh[0][1] = 0;
	gh[1][1] = 0;
	gh[2][1] = -fvr_re[0]*sig_xy;

	gh[3][1] = 0;
	gh[4][1] = 0;
	gh[5][1] = -fvr_re[1]*sig_xy;

	gh[6][1] = 0;
	gh[7][1] = 0;
	gh[8][1] = -fvr_re[2]*sig_xy;

	gh[9][1]  = 0;
	gh[10][1] = 0;
	gh[11][1] = -fvr_re[3]*sig_xy;
}

void LocalMatrixRZ::CalcLocalVectorHarmDiv(double *val_r_re, double *val_r_im)
{
	double fvr_re[4];
	double fvr_im[4];

	Mult_Plot((double*)c, val_r_re, fvr_re, 4);
	Mult_Plot((double*)c, val_r_im, fvr_im, 4);

	gh[0][0] = 0; 
	gh[1][0] = 0; 
	gh[2][0] = fvr_im[0]*sigma;

	gh[3][0] = 0;
	gh[4][0] = 0;
	gh[5][0] = fvr_im[1]*sigma;

	gh[6][0] = 0;
	gh[7][0] = 0;
	gh[8][0] = fvr_im[2]*sigma;

	gh[9][0]  = 0;
	gh[10][0] = 0;
	gh[11][0] = fvr_im[3]*sigma;

	gh[0][1] = 0;
	gh[1][1] = 0;
	gh[2][1] = -fvr_re[0]*sigma;

	gh[3][1] = 0;
	gh[4][1] = 0;
	gh[5][1] = -fvr_re[1]*sigma;

	gh[6][1] = 0;
	gh[7][1] = 0;
	gh[8][1] = -fvr_re[2]*sigma;

	gh[9][1]  = 0;
	gh[10][1] = 0;
	gh[11][1] = -fvr_re[3]*sigma;
}

void LocalMatrixRZ::CalcLocalVectorHarmTest(double *val_r_re, double *val_z_re, double *val_v_re,
										double *val_r_im, double *val_z_im, double *val_v_im)
{
	int i, j;

	double fr_re[4];
	double fz_re[4];
	double fv_re[4];

	double fr_im[4];
	double fz_im[4];
	double fv_im[4];

	Mult_Plot((double*)c, val_r_re, fr_re,  4);
	Mult_Plot((double*)c, val_z_re, fz_re,  4);
	Mult_Plot((double*)c, val_v_re, fv_re,  4);

	Mult_Plot((double*)c, val_r_im, fr_im,  4);
	Mult_Plot((double*)c, val_z_im, fz_im,  4);
	Mult_Plot((double*)c, val_v_im, fv_im,  4);

	for (i=0; i<4; i++)
	{
		gh[i*3][0]   = fr_re[i];
		gh[i*3+1][0] = fz_re[i];
		gh[i*3+2][0] = fv_re[i];

		gh[i*3][1]   = fr_im[i];
		gh[i*3+1][1] = fz_im[i];
		gh[i*3+2][1] = fv_im[i];
	}
}

void LocalMatrixRZ::CalcLocalVectorM1(double *val_r, double *val_z, double *val_v, double beta)
{
	int i, j;
	double f[12];

	for (int i=0; i<4; i++)
	{
		f[i*3]   = val_r[i];
		f[i*3+1] = val_z[i];
		f[i*3+2] = val_v[i];
	}

	Mult_Plot((double*)ac, f, g, 12);
}

void LocalMatrixRZ::CalcLocalVector(double *val)
{
	Mult_Plot((double*)c, val, g, 4);
}

int LocalMatrixRZ::GetStiffnessMatrix(const double &r0, const double &z0,const double &r3, const double &z3)
{
	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;
	double rp11, rp12, rp22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	rp11=(1+rk/h1)*(1+rk/h1)*log(1+h1/rk)-rk/h1-1.5; 
	rp12=-(1+rk/h1)*(rk/h1)*log(1+h1/rk)+rk/h1+0.5; 
	rp22=0.5-rk/h1+(rk/h1)*(rk/h1)*log(1+h1/rk);

	b_r2[2][2]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r11	+	rp11*z11;
	b_r2[2][3]=b_r2[3][2]=-(1./h1)*(rk+h1/2)*z11+(1./h2)*r12	+	rp12*z11;
	b_r2[2][0]=b_r2[0][2]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r11	+	rp11*z12;
	b_r2[2][1]=b_r2[1][2]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12	+	rp12*z12;

	b_r2[3][3]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r22	+	rp22*z11;
	b_r2[3][0]=b_r2[0][3]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12	+	rp12*z12;
	b_r2[3][1]=b_r2[1][3]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r22	+	rp22*z12;

	b_r2[0][0]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r11	+	rp11*z22;
	b_r2[0][1]=b_r2[1][0]=-(1./h1)*(rk+h1/2)*z22+(1./h2)*r12	+	rp12*z22;

	b_r2[1][1]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r22	+	rp22*z22;

	return 0;
}

int LocalMatrixRZ::GetMassMatrix(const double &r0, const double &z0, const double &r3, const double &z3)
{
	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	c[2][2]=		r11*z11;
	c[2][3]=c[3][2]=r12*z11;
	c[2][0]=c[0][2]=r11*z12;
	c[2][1]=c[1][2]=r12*z12;

	c[3][3]=		r22*z11;
	c[3][0]=c[0][3]=r12*z12;
	c[3][1]=c[1][3]=r22*z12;

	c[0][0]=		r11*z22;
	c[0][1]=c[1][0]=r12*z22;

	c[1][1]=		r22*z22;

	return 0; 
}

int LocalMatrixRZ::GetStiffnessMatrix0(const double &r0, const double &z0,const double &r3, const double &z3)
{
	double h1=r3-r0, h2=z0-z3, rk=r0;
	double r11, r12, r22, z11, z12, z22;

	r11=(h1/3)*(rk+h1/4); r12=(h1/6)*(rk+h1/2); r22=(h1/3)*(rk+3*h1/4);
	z11=h2/3; z12=h2/6; z22=h2/3;

	b[2][2]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r11;
	b[2][3]=b[3][2]=-(1./h1)*(rk+h1/2)*z11+(1./h2)*r12;
	b[2][0]=b[0][2]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r11;
	b[2][1]=b[1][2]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12;

	b[3][3]=		(1./h1)*(rk+h1/2)*z11+(1./h2)*r22;
	b[3][0]=b[0][3]=-(1./h1)*(rk+h1/2)*z12-(1./h2)*r12;
	b[3][1]=b[1][3]=(1./h1)*(rk+h1/2)*z12-(1./h2)*r22;

	b[0][0]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r11;
	b[0][1]=b[1][0]=-(1./h1)*(rk+h1/2)*z22+(1./h2)*r12;

	b[1][1]=		(1./h1)*(rk+h1/2)*z22+(1./h2)*r22;

	return 0;
}

int LocalMatrixRZ::GetStiffnessMatrix0(const double &r0, const double &z0,const double &r3, const double &z3,
									   const double &sig_xy, const double &sig_z)
{
	double t3 = 1/hr/hz;
	double t4 = hr*hr;
	double t6 = sig_z*t4*hr;
	double t8 = hz*hz;
	double t9 = sig_xy*rk*t8;
	double t10 = 4.0*t9;
	double t12 = sig_xy*hr*t8;
	double t13 = 2.0*t12;
	double t15 = sig_z*rk*t4;
	double t16 = 4.0*t15;
	double t19 = t3*(t6+t10+t13+t16)/12.0;
	double t20 = 2.0*t15;
	double t23 = t3*(t6-t10-t13+t20)/12.0;
	double t24 = 2.0*t9;
	double t27 = t3*(t6-t24-t12+t16)/12.0;
	double t30 = t3*(t6+t24+t12+t20)/12.0;
	double t31 = 3.0*t6;
	double t34 = t3*(t31+t10+t13+t16)/12.0;
	double t37 = t3*(t31-t24-t12+t16)/12.0;

	ba[0][0] = t19;
	ba[0][1] = t23;
	ba[0][2] = -t27;
	ba[0][3] = -t30;

	ba[1][0] = t23;
	ba[1][1] = t34;
	ba[1][2] = -t30;
	ba[1][3] = -t37;

	ba[2][0] = -t27;
	ba[2][1] = -t30;
	ba[2][2] = t19;
	ba[2][3] = t23;

	ba[3][0] = -t30;
	ba[3][1] = -t37;
	ba[3][2] = t23;
	ba[3][3] = t34;

	return 0;
}
