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
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin                                         
 *  Novosibirsk State Technical University,                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 December 10, 2024                                            
*/

#pragma once

class Rect_Local_Matrix
{
public:

	double J[2][2];
	double J_1_T[2][2];
	int type_of_rect;
	int alpha;
	double a[8][8];
	double b[4][4];
	double c[4][4];
	double g[8];
	double g_sin[4];
	double g_cos[4];
	double det_J;
	double det_J_abs;
	double x[4];
	double y[4];
	double a0_sin_y[4];
	double a0_cos_y[4];
	double h0_sin_y[4];
	double h0_cos_y[4];
	double *coords_1d; 
	double *sin_1d;
	double *cos_1d;
	long n_1d;
	double mu;
	double sigma;
	double sigma0;
	double omega;
	double sigma_omega;
	double mu_omega;
	Rect_Local_Matrix(long n_of_rect, double (*xy)[2], long (*nvtr)[4], long *type_of_rect,
		double *mu, double *sigma, double omega, long *nvkat, int alpha, double *sigma0,
		long n_1d, double *coords_1d, double *sin_1d, double *cos_1d);
	Rect_Local_Matrix(long n_of_rect, double (*xy)[2], long (*nvtr)[4]);
	Rect_Local_Matrix();
	~Rect_Local_Matrix();
	double l0(double x);
	double l1(double x);
	double Phi(long i, double x, double y);
	double D_phi(long i, long j, double xi, double eta);
	void Calc_J(int n_of_point);
	void Calc_J(double x, double y);
	void Calc_local_matrix();
	void Calc_elements_of_local_matrix_for_rectangle();
	void Calc_elements_of_local_matrix_for_quadrilateral(); 
	void Calc_Right_Part_From_1d();
};
