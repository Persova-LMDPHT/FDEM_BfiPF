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
 *  This file contains subroutines for basic vector and matrix-vector operations
 *                                                                                         
 *  Written by D.Sc. Denis V. Vagin, Prof. Marina G. Persova                       
 *  Novosibirsk State Technical University,                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)               
 * Version 2.0 December 10, 2024                                                                                                                      
 */                                                                                         

#include "stdafx.h"
#include "base_solver.h"
#include "ControlOMP.h"
extern ControlOMP omp;
extern ofstream logfile;
//-----------------------------------------------------------------
// Constuctor
//-----------------------------------------------------------------
Base_solver::Base_solver()
{
}
//-----------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------
Base_solver::~Base_solver()
{
}
//-----------------------------------------------------------------
// dot product
//-----------------------------------------------------------------
inline double Base_solver::Scal(double *x, double *y, int n)
{
	double sum = 0.0;
	int i;
	
	#pragma omp parallel shared(x, y) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for reduction(+:sum)
		for(i=0; i<n; i++)
		{
			sum += x[i]*y[i];
		}
	}

	return sum;
}
//-----------------------------------------------------------------
// vector euclidean norm
//-----------------------------------------------------------------
double Base_solver::Norm_Euclid(double *x, int n)
{
	double temp;

	temp = Scal(x,x,n);

	if(temp>0.0) return sqrt(temp);

	return 0.0;
}
//----------------------------------------------------------------------
// projection of a vector onto an axis
//----------------------------------------------------------------------
double Base_solver::Projection(double *vec, double *axis)
{
	return Scal(vec,axis,3)/Norm_Euclid(axis,3);
}
//---------------------------------------------------------------------------
// Matrix-vector multiplication  (in dense format)
//---------------------------------------------------------------------------
void Base_solver::Mult_Plot(double *a,double *pr,double *rez,int n)
{
	int i;
	for(i=0; i<n; i++)
		rez[i] = this->Scal(&a[i*n], pr, n);
}
//----------------------------------------------------------
// Givens rotations
//----------------------------------------------------------
int Base_solver::Givens1(double& x, double& y, double& c, double& s)
{
	double t = sqrt(x*x + y*y);

	if(fabs(t)<1e-30)
	{
		printf("Division by zero in Givens1.\n");
		return 1;
	}

	s = -y/t;
	c = x/t;
	x = x*c - y*s;
	y = 0.0;

	return 0;
}
//----------------------------------------------------------
void Base_solver::Givens2(double& x, double& y, double c, double s)
{
	double tempx = x;

	x = x*c - y*s;
	y = tempx*s + y*c;
}
//----------------------------------------------------------
int Base_solver::Givens(double *a, double *f, int n)
{
	int i, j;
	double c, s;

	for(i=0; i<n; i++)
	{
		if(Givens1(a[i*n+i], a[(i+1)*n+i], c, s)!=0)
			return 1;
		for(j=i+1; j<n; j++)
		{
			Givens2(a[i*n+j], a[(i+1)*n+j], c, s);
		}
		Givens2(f[i], f[i+1], c, s);
	}

	return 0;
}
//------------------------------------------------------------
// Solution of a system with a lower triangular matrix in a dense format
//------------------------------------------------------------
int Base_solver::Undirect(double *a, double *b, double *x, int n)
{
	int k, j;
	double s;

	for(k=n-1; k>=0; k--)
	{
		s = 0.0;
		for(j=k+1; j<n; j++)
		{
			s += a[k*n+j]*x[j];
		}
		x[k] = (b[k] - s)/a[k*n + k];
	}

	return 0;
}
//------------------------------------------------------------
// Solution of a SLAE with a square matrix whose lower triangle contains only one non - zero subdiagonal.	
// This is necessary if the Arnoldi orthogonalization fails
//------------------------------------------------------------
int Base_solver::Solve_square_subdiag(double *a, double *b, double *x, int n)
{
	int j, k;
	double mult;

	// reduce to upper triangular form
	for(j=0; j<n-1; j++)
	{
		mult = -a[(j+1)*n+j]/a[j*n+j];
		for(k=j+1; k<n; k++)
		{
			a[(j+1)*n+k] += mult*a[j*n+k];
		}
		b[j+1] += mult*b[j];
	}

	// a back
	Undirect(a, b, x, n);

	return 0;
}
//------------------------------------------------------------
// writing to the file the relative discrepancy with which they came out, eps, the number of iterations and the time of solving the SLAE
//------------------------------------------------------------
int Base_solver::WriteKitChrono(char *fname, double residual, double eps, int iter, double time)
{
	FILE *fp;
	printf("time=%lf\n", time);

	// kit
	if ((fp = fopen(fname, "w")) == 0)
	{
		printf("cannot open file %s\n", fname);
		return 1;
	}
	else
	{
		fprintf(fp, "r=%e\teps=%e\tkit=%d\t", residual, eps, iter);
		fprintf(fp, "time=%lf sec\n", time);
		fclose(fp);
	}

	return 0;
}
//------------------------------------------------------------                                                                                 
// writing to the file the relative discrepancy with which they came out, eps, the number of iterations and the time of solving the SLAE       
//------------------------------------------------------------                                                                                 
int Base_solver::Write_kit(char *fname, double residual, double eps, int iter, int time)
{
	int hours, minutes, seconds;
	FILE *fp;

	hours = time/3600;
	minutes = (time - hours*3600)/60;
	seconds = time - hours*3600 - minutes*60;

	printf("time=%ldh %ldm %lds %ldms\n", hours, minutes, seconds, time);

	// kit
	if((fp=fopen(fname,"w"))==0)
	{
		printf("cannot open file %s\n", fname);
		return 1;
	}
	else
	{
		fprintf(fp, "r=%e\teps=%e\tkit=%d\t", residual, eps, iter);
		fprintf(fp, "time=%ldh %ldm %lds %ldms\n", hours, minutes, seconds, time);
		fclose(fp);
	}	

	return 0;
}
//------------------------------------------------------------
int Base_solver::Write_kit(char *fname, double residual, double eps, int iter, int time, int n)
{
	int hours, minutes, seconds;
	FILE *fp;

	hours = time/3600;
	minutes = (time - hours*3600)/60;
	seconds = time - hours*3600 - minutes*60;

	printf("time=%ldh %ldm %lds\n", hours, minutes, seconds);

	// kit
	if((fp=fopen(fname, "a"))==0)
	{
		printf("cannot open file %s\n", fname);
		return 1;
	}
	else
	{
		fprintf(fp, "r=%e\teps=%e\tkit=%d\t", residual, eps, iter);
		fprintf(fp, "time=%ldh %ldm %lds\t", hours, minutes, seconds);
		fprintf(fp, "n=%d\n", n);
		fclose(fp);
	}	

	return 0;
}
//------------------------------------------------------------
int Base_solver::Write_kit(char *fname, double residual, double eps, int iter, int time, int n, int size_jg)
{
	int hours, minutes, seconds;
	FILE *fp;

	hours = time/3600;
	minutes = (time - hours*3600)/60;
	seconds = time - hours*3600 - minutes*60;

	printf("time=%ldh %ldm %lds\n", hours, minutes, seconds);

	// kit
	if((fp=fopen(fname, "a"))==0)
	{
		printf("cannot open file %s\n", fname);
		return 1;
	}
	else
	{
		fprintf(fp, "r=%e\teps=%e\tkit=%d\t", residual, eps, iter);
		fprintf(fp, "time=%ldh %ldm %lds\t", hours, minutes, seconds);
		fprintf(fp, "n=%d\t", n);
		fprintf(fp, "size_jg=%d\n", size_jg);
		fclose(fp);
	}	

	return 0;
}
//------------------------------------------------------------
int Base_solver::Write_kit(char *fname, double residual, double eps, int iter, int time, double change_of_solution)
{
	int hours, minutes, seconds;
	FILE *fp;

	hours = time/3600;
	minutes = (time - hours*3600)/60;
	seconds = time - hours*3600 - minutes*60;

	printf("time=%ldh %ldm %lds\n", hours, minutes, seconds);

	// kit
	if((fp=fopen(fname, "a"))==0)
	{
		printf("cannot open file %s\n", fname);
		return 1;
	}
	else
	{
		fprintf(fp, "r=%e\teps=%e\tkit=%d\t", residual, eps, iter);
		fprintf(fp, "time=%ldh %ldm %lds\t", hours, minutes, seconds);
		fprintf(fp, "change of solution=%e\n", change_of_solution);
		fclose(fp);
	}	

	return 0;
}
//------------------------------------------------------------
// relative error
//------------------------------------------------------------
double Base_solver::Relative_Error(double *analytic, double *numeric, int n)
{
	double *razn;
	double norm_analytic;
	double norm_razn;
	double error;
	int i;

	razn = new double[n];
	if(razn == 0)
	{
		char var[] = {"razn"};
		char func[] = {"Relative_Error"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	for(i=0;i<n;i++)
		razn[i] = analytic[i] - numeric[i];

	norm_analytic = Norm_Euclid(analytic, n);
	norm_razn = Norm_Euclid(razn, n);

	error = norm_razn/norm_analytic;

	delete [] razn;

	return error;
}
//-----------------------------------------------------------
// linear interpolation
//-----------------------------------------------------------
double Base_solver::Spline(double x, int n, double *xyz, double *values)
{
	double s, xi;
	int i, t, flag;

	flag = 0;

	for(i=0; i<n; i++)
	{
		if(x >= xyz[i]  &&  x <= xyz[i+1])
		{
			t = i;
			flag = 1;
			break;
		}
	}

		if(flag == 1)
		{
			xi = (x - xyz[t])/(xyz[t+1] - xyz[t]);
			s = (1.0 - xi)*values[t] + xi*values[t+1];
		}
		else
		{
			s = 0.0;
		}

		return s;
}
//------------------------------------------------------------------------
// complex conjugate dot product
//------------------------------------------------------------------------
std::complex<double> Base_solver::ScalCmplx(double *x, double *y, int nb)
{
	double s_re = 0;
	double s_im = 0;
	int i;

	#pragma omp parallel shared(x, y) private(i) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for reduction(+:s_re) reduction(+:s_im)
		for (i=0; i<nb; i++)
		{
			s_re += x[i*2]*y[i*2]   - x[i*2+1]*y[i*2+1];
			s_im += x[i*2]*y[i*2+1] + x[i*2+1]*y[i*2];
		}
	}

	return std::complex<double>(s_re, s_im);
}
//------------------------------------------------------------------------
// dot product for complex-valued vectors
//------------------------------------------------------------------------
std::complex<double> Base_solver::ScalCmplxTrue(double *x, double *y, int nb)
{
	double s_re = 0;
	double s_im = 0;
	int i;

	for (i=0; i<nb; i++)
	{
		s_re += x[i*2]*y[i*2]   - x[i*2+1]*y[i*2+1];
		s_im += x[i*2]*y[i*2+1] + x[i*2+1]*y[i*2];
	}

	return std::complex<double>(s_re, s_im);
}
//------------------------------------------------------------------------
// multiplication of a vector by a complex number
//------------------------------------------------------------------------
void Base_solver::MultCmplxNumVect(std::complex<double> a, double *x, double *y, int nb)
{
	int i;

	for (i=0; i<nb; i++)
	{
		y[i*2] = x[i*2]*real(a) - x[i*2+1]*imag(a);
		y[i*2+1] = x[i*2]*imag(a) + x[i*2+1]*real(a);
	}
}
//------------------------------------------------------------------------
// z = x + a*y
//------------------------------------------------------------------------
void Base_solver::Cmplx_axpy(std::complex<double> a, double *x, double *y, double *z, int nb)
{
	int i;
	double aRe = real(a);
	double aIm = imag(a);
	double t_re, t_im;

	#pragma omp parallel shared(x, y, z) private(i, t_re, t_im) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for 
		for (i=0; i<nb; i++)
		{
			t_re = x[i*2]   + (y[i*2]*aRe - y[i*2+1]*aIm);
			t_im = x[i*2+1] + (y[i*2]*aIm + y[i*2+1]*aRe);

			z[i*2]   = t_re;
			z[i*2+1] = t_im;
		}
	}
}
//------------------------------------------------------------------------
// multiplication of components of one complex vector by components of another complex vector
//------------------------------------------------------------------------
void Base_solver::MultCmplxVectVect(int nb, double *a, double *b, double *c)
{
	int i;
	std::complex<double> x, y, z;

	#pragma omp parallel shared(a, b, c) private(x, y, z) num_threads(omp.GetNumberOfThreads()) 
	{
		#pragma omp for 
		for (i=0; i<nb; i++)
		{
			x = std::complex<double>(a[i*2], a[i*2+1]);
			y = std::complex<double>(b[i*2], b[i*2+1]);
			z = x*y;
			c[i*2]   = real(z);
			c[i*2+1] = imag(z);
		}
	}
}
//------------------------------------------------------------------------
// division of the components of one complex vector into components of another complex vector
//------------------------------------------------------------------------
void Base_solver::DivCmplxVectVect(int nb, double *a, double *b, double *c)
{
	int i;
	std::complex<double> x, y, z;

	for (i=0; i<nb; i++)
	{
		x = std::complex<double>(a[i*2], a[i*2+1]);
		y = std::complex<double>(b[i*2], b[i*2+1]);
		z = x/y;
		c[i*2]   = real(z);
		c[i*2+1] = imag(z);
	}
}
//------------------------------------------------------------------------
double Base_solver::GetMinimiz(int n, double *ap, double *r)
{
	double w = 0;
	double rAp = Scal(ap, r, n);

	if (fabs(rAp) > 0)
		w = Scal(ap, r, n)/rAp;

	return w;
}
