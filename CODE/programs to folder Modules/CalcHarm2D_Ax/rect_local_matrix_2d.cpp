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
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024
*/

#include "stdafx.h"
#include "rect_local_matrix_2d.h"
#include "gauss3_2d.h"
#include "For_Solvers.h"

Rect_Local_Matrix::Rect_Local_Matrix()
{
}

Rect_Local_Matrix::Rect_Local_Matrix(long n_of_rect, double (*xy)[2], long (*nvtr)[4])
{
	long i;

	for(i=0; i<4; i++)
	{
		this->x[i] = xy[nvtr[n_of_rect][i]][0];
		this->y[i] = xy[nvtr[n_of_rect][i]][1];
	}

	for(i=0; i<8; i++)
		g[i] = 0.0;
}

Rect_Local_Matrix::Rect_Local_Matrix(long n_of_rect, double (*xy)[2], long (*nvtr)[4],
									 long *type_of_rect, double *mu, double *sigma,
									 double omega, long *nvkat, int alpha, double *sigma0,
									 long n_1d, double *coords_1d, double *sin_1d,
									 double *cos_1d)
{
	long i;

	this->alpha = alpha;

	this->mu = mu[nvkat[n_of_rect]];
	this->sigma = sigma[nvkat[n_of_rect]];
	this->sigma0 = sigma0[nvkat[n_of_rect]];
	this->omega = omega;

	this->n_1d = n_1d;
	this->coords_1d = coords_1d;
	this->sin_1d = sin_1d;
	this->cos_1d = cos_1d;

	for(i=0; i<4; i++)
	{
		this->x[i] = xy[nvtr[n_of_rect][i]][0];
		this->y[i] = xy[nvtr[n_of_rect][i]][1];
	}

	this->type_of_rect = type_of_rect[n_of_rect];
}

Rect_Local_Matrix::~Rect_Local_Matrix()
{
}

void Rect_Local_Matrix::Calc_J(int n_of_point)
{
	long i;

	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

	for(i=0; i<4; i++)
	{
		J[0][0] += x[i]*gauss_3_d_phi_rect[n_of_point][i][0];
		J[0][1] += x[i]*gauss_3_d_phi_rect[n_of_point][i][1];
		J[1][0] += y[i]*gauss_3_d_phi_rect[n_of_point][i][0];
		J[1][1] += y[i]*gauss_3_d_phi_rect[n_of_point][i][1];
	}

	this->det_J = J[0][0]*J[1][1] - J[1][0]*J[0][1];
	this->det_J_abs = fabs(det_J);

	J_1_T[0][0] = J[1][1]/det_J;
	J_1_T[1][0] = -J[0][1]/det_J;
	J_1_T[0][1] = -J[1][0]/det_J;
	J_1_T[1][1] = J[0][0]/det_J;
}
//---------------------------------------------------
void Rect_Local_Matrix::Calc_J(double x0, double y0)
{
	long i;

	J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;

	for(i=0; i<4; i++)
	{
		J[0][0] += x[i]*D_phi(i, 0, x0, y0);
		J[0][1] += x[i]*D_phi(i, 1, x0, y0);
		J[1][0] += y[i]*D_phi(i, 0, x0, y0);
		J[1][1] += y[i]*D_phi(i, 1, x0, y0);
	}

	this->det_J = J[0][0]*J[1][1] - J[1][0]*J[0][1];
	this->det_J_abs = fabs(det_J);

	J_1_T[0][0] = J[1][1]/det_J;
	J_1_T[1][0] = -J[0][1]/det_J;
	J_1_T[0][1] = -J[1][0]/det_J;
	J_1_T[1][1] = J[0][0]/det_J;
}

void Rect_Local_Matrix::Calc_elements_of_local_matrix_for_rectangle()
{

	double hx, hy, hy_hx, hx_hy, hy_hx_2, hx_hy_2;
	double t;

	hx = x[1] - x[2];
	hy = y[1] - y[2];

	hx_hy = hx/hy/3.0;
	hx_hy_2 = hx_hy*0.5;
	hy_hx = hy/hx/3.0;
	hy_hx_2 = hy_hx*0.5;

	b[0][0] = b[1][1] = b[2][2] = b[3][3] =  hy_hx + hx_hy;
	b[0][1] = b[2][3] = b[1][0] = b[3][2] = -hy_hx + hx_hy_2;
	b[0][2] = b[1][3] = b[2][0] = b[3][1] = hy_hx_2 - hx_hy;
	b[0][3] = b[1][2] = b[3][0] = b[2][1] = -hy_hx_2 - hx_hy_2;	

	t = hx*hy/36.0;

	c[0][0] = c[1][1] = c[2][2] = c[3][3] = 4.0*t;

	c[0][1] = c[0][2] = c[1][3] = c[2][3] =
	c[1][0] = c[2][0] = c[3][1] = c[3][2] = 2.0*t;
	c[0][3] = c[1][2] = c[2][1] = c[3][0] = t;
}

void Rect_Local_Matrix::Calc_elements_of_local_matrix_for_quadrilateral()
{

	int i, j, i1, j1;
	double gauss_3_mult;
	double grad_all[4][2];

	for(i=0; i<4; i++)
	for(j=0; j<4; j++)
		b[i][j] = c[i][j] = 0.0;

	for(i=0; i<9; i++)
	{
		this->Calc_J(i);
		gauss_3_mult = gauss_3_A_all_rect[i]*det_J_abs;

		for(j=0; j<4; j++)
			Mult_Plot((double*)J_1_T, (double*)gauss_3_d_phi_rect[i][j], (double*)grad_all[j], 2);			

		for(i1=0; i1<4; i1++)
		for(j1=0; j1<4; j1++)
		{
			c[i1][j1] += gauss_3_phi_rect[i][i1]*gauss_3_phi_rect[i][j1]*gauss_3_mult;
			b[i1][j1] += Scal((double*)grad_all[i1],(double*)grad_all[j1], 2)*gauss_3_mult;
		}
	}    
}

void Rect_Local_Matrix::Calc_local_matrix()
{
	int i, j;
	double tmp;
	double sigma0_sigma_omega;
	double sigma0_1_sigma;

	Calc_Right_Part_From_1d();
	for(i=0; i<8; i++)
		g[i] = 0.0;

	if(type_of_rect < 5)
	{
		Calc_elements_of_local_matrix_for_rectangle();
	}
	else
	{
		Calc_elements_of_local_matrix_for_quadrilateral();
	}

	if(alpha == 0)
	{
		sigma_omega = sigma*omega;

		for(i=0; i<4; i++)
		for(j=0; j<4; j++)
		{
			tmp = sigma_omega*c[i][j];
			a[i*2][j*2] = a[i*2+1][j*2+1] = b[i][j]/mu;
			a[i*2][j*2+1] = -tmp;
			a[i*2+1][j*2] = tmp;
		}

		Mult_Plot((double*)c, a0_cos_y, g_sin, 4);
		Mult_Plot((double*)c, a0_sin_y, g_cos, 4);

		sigma0_sigma_omega = (sigma0 - sigma)*omega;

		for(i=0; i<4; i++)
		{
			g[i*2]   = -g_sin[i]*sigma0_sigma_omega;
			g[i*2+1] = g_cos[i]*sigma0_sigma_omega;
		}
	}
	else
	{
		mu_omega = mu*omega;

		for(i=0; i<4; i++)
		for(j=0; j<4; j++)
		{
			tmp = mu_omega*c[i][j];
			a[i*2][j*2] = a[i*2+1][j*2+1] = b[i][j]/sigma;
			a[i*2][j*2+1] = -tmp;
			a[i*2+1][j*2] = tmp;
		}

		Mult_Plot((double*)b, h0_sin_y, g_sin, 4);
		Mult_Plot((double*)b, h0_cos_y, g_cos, 4);

		sigma0_1_sigma = (1.0/sigma0  - 1.0/sigma);

		for(i=0; i<4; i++)
		{
			g[i*2]   = g_sin[i]*sigma0_1_sigma;
			g[i*2+1] = g_cos[i]*sigma0_1_sigma;
		}
	}
}

double Rect_Local_Matrix::l0(double x)
{
	return (1.0 - x)*0.5;
}

double Rect_Local_Matrix::l1(double x)
{
	return (x + 1.0)*0.5;
}

double Rect_Local_Matrix::Phi(long i, double x, double y)
{
	switch(i)
	{
	case 0:
		return l0(x)*l1(y);
		break;
	case 1:
		return l1(x)*l1(y);
		break;
	case 2:
		return l0(x)*l0(y);
		break;
	default:
		return l1(x)*l0(y);
	}
}

double Rect_Local_Matrix::D_phi(long i, long j, double x, double y)
{
	if(j==0)
	{
		switch(i)
		{
			case 0:  return -y/4.0-1.0/4.0;
			case 1:  return	y/4.0+1.0/4.0;
			case 2:  return -1.0/4.0+y/4.0;
			default: return	1.0/4.0-y/4.0;
		}
	}
	else
	{
		switch(i)
		{
			case 0:  return 1.0/4.0-x/4.0;
			case 1:  return x/4.0+1.0/4.0;
			case 2:  return -1.0/4.0+x/4.0;
			default: return -x/4.0-1.0/4.0;
		}
	}
}

void Rect_Local_Matrix::Calc_Right_Part_From_1d()
{
	int i;

	if(alpha == 0)
	{
		for(i=0; i<4; i++)
		{
			a0_sin_y[i] = Spline(y[i], n_1d, coords_1d, sin_1d); 
			a0_cos_y[i] = Spline(y[i], n_1d, coords_1d, cos_1d); 
		}
	}
	else
	{
		for(i=0; i<4; i++)
		{
			h0_sin_y[i] = Spline(y[i], n_1d, coords_1d, sin_1d); 
			h0_cos_y[i] = Spline(y[i], n_1d, coords_1d, cos_1d); 
		}
	}
}
