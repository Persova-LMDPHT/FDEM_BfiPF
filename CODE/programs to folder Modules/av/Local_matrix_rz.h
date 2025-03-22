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
 *  This file contains headers of functions for calculating 2D local matirxes
 * 
 * 
 *  Written by  Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 * Version 2.0 December 10, 2024                                            
*/

#pragma once

class LocalMatrixRZ
{
public:
	double b[4][4];
	double ba[4][4];
	double b_r2[4][4];
	double c[4][4];
	double dr[4][4];
	double dz[4][4];

	double ab[12][12];
	double ac[12][12];
	double aa[12][12];

	double ah[12][12][2];
	double gh[12][2];

	double f[4];
	double g[12];

	double rk;
	double zk;
	double hr;
	double hz;

	double r0;
	double r1;
	double z0;
	double z1;

	double mu;
	double sigma;
	double beta;

	// Calculate the local stiffness matrix of a finite element
	void CalcLocalMatrixB();
	// Calculate the local stiffness matrix of a finite element with a term from the rotor
	void CalcLocalMatrixB_r2();
	// Calculate the local stiffness matrix of an anisotropic finite element
	void CalcLocalMatrixB(double sig_xy, double sig_z);
	// Calculate the local mass matrix of the finite element
	void CalcLocalMatrixC();

	// Calculate the local matrix from integrals with Psi*d_Psi/d_r of the finite element
	void CalcLocalMatrixDr();
	// Calculate the local matrix from the integrals with Psi*d_Psi/d_z of the finite element
	void CalcLocalMatrixDz();
	void CalcLocalMatrix();
	void CalcLocalMatrix2(double beta);
	void CalcLocalVector(double *val_r, double *val_z, double *val_v, double beta);
	void CalcLocalVectorM(double *val_r, double *val_z, double *val_v, double beta);
	void CalcLocalVectorM1(double *val_r, double *val_z, double *val_v, double beta);
	void CalcLocalVector(double *val);

	void CalcLocalMatrixHarm(double omega);
	void CalcLocalMatrixHarm(double omega, double sig_xy, double sig_z);
	void CalcLocalVectorHarm(double *val_r_re, double *val_z_re, double *val_v_re,
		double *val_r_im, double *val_z_im, double *val_v_im);

	void CalcLocalVectorHarmTest(double *val_r_re, double *val_z_re, double *val_v_re,
		double *val_r_im, double *val_z_im, double *val_v_im);

	void CalcLocalVectorHarmDiv(double *val_r_re, double *val_r_im, double sig_xy, double sig_z);
	void CalcLocalVectorHarmDiv(double *val_r_re, double *val_r_im);

	int GetMassMatrix(const double &r0, const double &z0,const double &r3, const double &z3);
	int GetStiffnessMatrix0(const double &r0, const double &z0,const double &r3, const double &z3);
	int GetStiffnessMatrix0(const double &r0, const double &z0,const double &r3, const double &z3,const double &sig_xy, const double &sig_z);
	int GetStiffnessMatrix(const double &r0, const double &z0,const double &r3, const double &z3);

	LocalMatrixRZ(double r0, double r1, double z0, double z1, double mu, double sigma, double beta);
};
