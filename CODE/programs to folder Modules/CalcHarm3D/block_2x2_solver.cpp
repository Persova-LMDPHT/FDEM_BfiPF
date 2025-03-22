/*
 * GENERAL REMARKS                                                                           
 *                                                                                           
 *  This code is freely available under the following conditions:                            
 *                                                                                           
 *  1) The code is to be used only for non-commercial purposes.                              
 *  2) No changes and modifications to the code without prior permission of the developer.   
 *  3) No forwarding the code to a third party without prior permission of the developer.    
 *                                                                                           
 *                       
 *  This file contains subroutines for matrix-vector operations with a matrix stored in 2x2-block sparse format          
 *                                                                                           
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/block_2x2_solver.cpp                                                        
 *  Novosibirsk State Technical University,                                                  
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                        
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                       
 *                                                                 
*/                                                                                           

#include "stdafx.h"
#include "block_2x2_solver.h"
#include "ControlOMP.h"
extern ControlOMP omp;
extern ofstream logfile;
//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------
Block_2x2_solver::Block_2x2_solver()
{
}
//-------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------
Block_2x2_solver::~Block_2x2_solver()
{
}
//------------------------------------------------------------------------
//  multiplication of a single block
//------------------------------------------------------------------------
inline void Block_2x2_solver::Mult_block_2x2(double *a, int size, double *x, double *y)
{
	if (size == 2)
	{
		y[0] += a[0]*x[0] - a[1]*x[1];
		y[1] += a[1]*x[0] + a[0]*x[1];
	} 
	else
	{
		y[0] += a[0]*x[0];
		y[1] += a[0]*x[1];
	}
}
//------------------------------------------------------------------------
// multiplication of a transposed block
//------------------------------------------------------------------------
inline void Block_2x2_solver::Mult_MV_block_2x2_transp(double *a, int size, double *x, double *y)
{
	if (size == 2)
	{
		y[0] +=  a[0]*x[0] + a[1]*x[1];
		y[1] += -a[1]*x[0] + a[0]*x[1];
	} 
	else
	{
		y[0] += a[0]*x[0];
		y[1] += a[0]*x[1];
	}
}
//------------------------------------------------------------------------
// multiplication of a block matrix by a vector
//------------------------------------------------------------------------
void Block_2x2_solver::Mult_MV_block_2x2(int nb, int *ig, int *jg, int *idi, int *ijg,
										 double *di_block, double *ggl_block, double *x, double *y, double *y_omp)
{
	int i, j, k;
	int size; // block size
	int ib, jb;
	
	int adr;
	int rank; // number of the current thread
	int threadNum = omp.GetNumberOfThreads();
	if (threadNum == 1 || nb < omp.GetNMinSparseMultMV2x2()) 
	{
		for(i=0; i<nb*2; i++) // vector-result is filled by 0
			y[i] = 0.0;

		for(i=0; i<nb; i++)
		{
			// diagonal blocks multiplication
			ib = i*2;
			size = idi[i+1] - idi[i];
			Mult_block_2x2(&di_block[idi[i]], size, &x[ib], &y[ib]);

			// non-diagonal blocks multiplication
			for(j=ig[i]; j<=ig[i+1]-1; j++)
			{
				jb = jg[j]*2;
				k = ijg[j];
				size = ijg[j+1] - k;
				Mult_block_2x2(&ggl_block[k], size, &x[jb], &y[ib]);
				Mult_block_2x2(&ggl_block[k], size, &x[ib], &y[jb]);
			}
		}
	}
	else
	{
		#pragma omp parallel shared(ig, jg, idi, ijg, di_block, ggl_block, x, y, y_omp) private(i, j, k, size, ib, jb, rank, adr) num_threads(threadNum)
		{
			#pragma omp for nowait
			for(i=0; i<nb*2; i++) // vector-result is filled by 0
			{
				y[i] = 0.0;
			}

			#pragma omp for
			for (i=0; i<nb*2*(threadNum-1); i++)
			{
				y_omp[i] = 0.0;
			}

			#pragma omp for
			for(i=0; i<nb; i++)
			{
				rank = omp_get_thread_num();
				if (rank == 0)
				{
					// diagonal blocks multiplication
					ib = i*2;
					size = idi[i+1] - idi[i];
					Mult_block_2x2(&di_block[idi[i]], size, &x[ib], &y[ib]);

					// non-diagonal blocks multiplication
					for(j=ig[i]; j<=ig[i+1]-1; j++)
					{
						jb = jg[j]*2;
						k = ijg[j];
						size = ijg[j+1] - k;
						Mult_block_2x2(&ggl_block[k], size, &x[jb], &y[ib]);
						Mult_block_2x2(&ggl_block[k], size, &x[ib], &y[jb]);
					}
				} 
				else
				{
					adr = (rank - 1)*nb*2;

					// diagonal blocks multiplication
					ib = i*2;
					size = idi[i+1] - idi[i];
					Mult_block_2x2(&di_block[idi[i]], size, &x[ib], &y_omp[adr + ib]);

					// non-diagonal blocks multiplication
					for(j=ig[i]; j<=ig[i+1]-1; j++)
					{
						jb = jg[j]*2;
						k = ijg[j];
						size = ijg[j+1] - k;
						Mult_block_2x2(&ggl_block[k], size, &x[jb], &y_omp[adr + ib]);
						Mult_block_2x2(&ggl_block[k], size, &x[ib], &y_omp[adr + jb]);
					}
				}// else
			}// i

			//reduction y_omp
			#pragma omp for
			for (i=0; i<nb*2; i++)
			{
				for (j=0; j<threadNum-1; j++)
					y[i] += y_omp[j*nb*2 + i];
			}	
		}// parallel
	}// else
}
//------------------------------------------------------------------------
// multiplication of a transposed block matrix by a vector
//------------------------------------------------------------------------
void Block_2x2_solver::Mult_MV_block_2x2_transp(int nb, int *ig, int *jg, int *idi, int *ijg,
										 double *di_block, double *ggl_block, double *x, double *y, double *y_omp)
{
	int i, j, k;
	int size; // block size
	int ib, jb;

	int adr;
	int rank; // number of the current thread
	int threadNum = omp.GetNumberOfThreads();

	if (threadNum == 1 || nb < omp.GetNMinSparseMultMV2x2()) 
	{
		for(i=0; i<nb*2; i++) // vector-result is filled by 0
			y[i] = 0.0;

		for(i=0; i<nb; i++)
		{
			// diagonal blocks multiplication
			ib = i*2;
			size = idi[i+1] - idi[i];
			Mult_MV_block_2x2_transp(&di_block[idi[i]], size, &x[ib], &y[ib]);

			// non-diagonal blocks multiplication
			for(j=ig[i]; j<=ig[i+1]-1; j++)
			{
				jb = jg[j]*2;
				k = ijg[j];
				size = ijg[j+1] - k;
				Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[jb], &y[ib]);
				Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[ib], &y[jb]);
			}
		}
	}
	else
	{
#pragma omp parallel shared(ig, jg, idi, ijg, di_block, ggl_block, x, y, y_omp) private(i, j, k, size, ib, jb, rank, adr) num_threads(threadNum)
		{
#pragma omp for nowait
			for(i=0; i<nb*2; i++) // vector-result is filled by 0
			{
				y[i] = 0.0;
			}

#pragma omp for
			for (i=0; i<nb*2*(threadNum-1); i++)
			{
				y_omp[i] = 0.0;
			}

#pragma omp for
			for(i=0; i<nb; i++)
			{
				rank = omp_get_thread_num();

				if (rank == 0)
				{
					// diagonal blocks multiplication
					ib = i*2;
					size = idi[i+1] - idi[i];
					Mult_MV_block_2x2_transp(&di_block[idi[i]], size, &x[ib], &y[ib]);

					// non-diagonal blocks multiplication
					for(j=ig[i]; j<=ig[i+1]-1; j++)
					{
						jb = jg[j]*2;
						k = ijg[j];
						size = ijg[j+1] - k;
						Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[jb], &y[ib]);
						Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[ib], &y[jb]);
					}
				} 
				else
				{
					adr = (rank - 1)*nb*2;

					// diagonal blocks multiplication
					ib = i*2;
					size = idi[i+1] - idi[i];
					Mult_MV_block_2x2_transp(&di_block[idi[i]], size, &x[ib], &y_omp[adr + ib]);

					// non-diagonal blocks multiplication
					for(j=ig[i]; j<=ig[i+1]-1; j++)
					{
						jb = jg[j]*2;
						k = ijg[j];
						size = ijg[j+1] - k;
						Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[jb], &y_omp[adr + ib]);
						Mult_MV_block_2x2_transp(&ggl_block[k], size, &x[ib], &y_omp[adr + jb]);
					}
				}// else
			}// i

			//reduction y_omp
#pragma omp for
			for (i=0; i<nb*2; i++)
			{
				for (j=0; j<threadNum-1; j++)
					y[i] += y_omp[j*nb*2 + i];
			}	
		}// parallel
	}// else
}
//------------------------------------------------------------------------
// Build block-diagonal preconditioner
//------------------------------------------------------------------------
int Block_2x2_solver::Build_block_diag_preconditioner(int nb, int *idi, 
	double *di_block, double *df, int *idi_f, double *ggl_f, double *ggu_f)
{
	int i;
	int adr; // the address of the beginning of storage of off-diagonal elements of the source block in di_block[]
	int adr_f; // the address of the beginning of storage of off-diagonal elements of the factorized block in ggl_f[], ggu_f[]
	double b, c;

	for(i=0; i<nb; i++)
	{
		adr = idi[i];

		if(idi[i+1] - adr == 2) // block size==2
		{
			b = di_block[adr];
			c = di_block[adr+1];

			// df
			df[i*2] = sqrt(b);
			df[i*2+1]=sqrt(b + c*c/b);

			// ggl_f, ggu_f
			adr_f = idi_f[i];
			ggl_f[adr_f] = c/sqrt(b);
			ggu_f[adr_f] = -ggl_f[adr_f];
		}
		else
		{//block size==1
			df[i*2]=df[i*2+1] = sqrt(b);
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Build block-diagonal preconditioner
//------------------------------------------------------------------------
int Block_2x2_solver::Build_block_diag_preconditioner(int nb, int *idi, 
													  double *di_block, double *df, double *ggl_f, double *ggu_f)
{
	int i;
	int adr; 
	double b, c;
	
	#pragma omp parallel shared(idi, di_block, df, ggl_f, ggu_f) private(i, adr, b, c) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<nb; i++)
		{
			adr = idi[i];

			if(idi[i+1] - adr == 2) // block size==2
			{

				b = di_block[adr];
				c = di_block[adr+1];

				df[i*2] = sqrt(b);

				ggl_f[i] = c/sqrt(b);
				ggu_f[i] = -c/sqrt(b);
				df[i*2+1]=sqrt(b - ggl_f[i]*ggu_f[i]);
			}
			else // block size==1
			{
				b = di_block[adr];
				df[i*2]=df[i*2+1] = sqrt(b);

				ggl_f[i] = ggu_f[i] = 0.0;
			}
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Build block-diagonal preconditioner
//------------------------------------------------------------------------
int Block_2x2_solver::Build_complex_diag_preconditioner(int nb, int *idi, double *di_block, double *df)
{
	int i;
	int adr; 
	std::complex<double> a, d;

	for(i=0; i<nb; i++)
	{
		adr = idi[i];

		if(idi[i+1] - adr == 2) //block size==2
		{
			a = std::complex<double>(di_block[adr], di_block[adr+1]);
			d = std::complex<double>(1, 0)/a;

			df[i*2]   = real(d);
			df[i*2+1] = imag(d);
		}
		else //block size==1
		{
			df[i*2]   = 1.0/di_block[adr];
			df[i*2+1] = 0.0;
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Solve lower triangular system with block-diagonal matrix
//------------------------------------------------------------------------
int Block_2x2_solver::solve_l_blockdiag(int nb, double *df, int *idi_f, double *ggl_f,
										double *f, double *x)
{
	int i, ib, ib1, size;

	for(i=0; i<nb; i++)
	{
		size = idi_f[i+1] - idi_f[i];

		ib = i*2;
		ib1 = i*2+1;

		if (size == 1)
		{
			x[ib] = f[ib]/df[ib];
			x[ib1] = (f[ib1] - x[ib1]*ggl_f[idi_f[i]])/df[ib1];
		} 
		else
		{
			x[ib] = f[ib]/df[ib];
			x[ib1] = f[ib1]/df[ib1];
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Solve upper triangular system with block-diagonal matrix
//------------------------------------------------------------------------
int Block_2x2_solver::solve_u_blockdiag(int nb, double *df, int *idi_f, double *ggu_f,
										double *f, double *x)
{
	int i, ib, ib1, size;

	for(i=0; i<nb; i++)
	{
		size = idi_f[i+1] - idi_f[i];

		ib = i*2;
		ib1 = i*2+1;

		if (size == 1)
		{
			x[ib1] = f[ib1]/df[ib1];
			x[ib] = (f[ib] - x[ib1]*ggu_f[idi_f[i]])/df[ib];
		} 
		else
		{
			x[ib1] = f[ib1]/df[ib1];
			x[ib] = f[ib]/df[ib];
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Solve lower triangular system with block-diagonal matrix
//------------------------------------------------------------------------
int Block_2x2_solver::solve_l_blockdiag(int nb, double *df, double *ggl_f, double *f, double *x)
{
	int i, ib, ib1;

	#pragma omp parallel shared(df, ggl_f, f, x) private(i, ib, ib1) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<nb; i++)
		{
			ib = i*2;
			ib1 = i*2+1;

			x[ib] = f[ib]/df[ib];
			x[ib1] = (f[ib1] - x[ib]*ggl_f[i])/df[ib1];
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Solve upper triangular system with block-diagonal matrix
//------------------------------------------------------------------------
int Block_2x2_solver::solve_u_blockdiag(int nb, double *df, double *ggu_f, double *f, double *x)
{
	int i, ib, ib1;

	#pragma omp parallel shared(df, ggu_f, f, x) private(i, ib, ib1) num_threads(omp.GetNumberOfThreads())
	{
		#pragma omp for
		for(i=0; i<nb; i++)
		{
			ib = i*2;
			ib1 = i*2+1;

			x[ib1] = f[ib1]/df[ib1];
			x[ib] = (f[ib] - x[ib1]*ggu_f[i])/df[ib];
		}
	}
	return 0;
}
//------------------------------------------------------------------------
// Complex-valued incomplete Cholesky factorization
//------------------------------------------------------------------------
int Block_2x2_solver::LLT_Cmplx(int nb, int *ig, int *jg, int *idi,
								int *ijg, double *di, double *gg,
								double *d, double *sg)
{
	int i, l, k, k1;
	std::complex<double> s, temp, a, b;
	int n = nb*2;
	int sz_sg = ig[nb]*2;
	int sz;

	cout << "Complex LLT begin...\n";
	logfile << "Complex LLT begin...\n";

	for(i=0; i<n; i++)
		d[i] = 0;

	for (i=0; i<sz_sg; i++)
		sg[i] = 0;

	for (i=0; i<nb; i++)
	{
		for(l=ig[i]; l<=ig[i+1]-1; l++)
		{
			s = std::complex<double>(0, 0);

			for(k=ig[i]; k<=l-1; k++)
			{
				for(k1=ig[jg[l]]; k1<=ig[jg[l]+1]-1; k1++)
				{
					if(jg[k1] == jg[k])
					{
						a = std::complex<double>(sg[k1*2], sg[k1*2+1]);
						b = std::complex<double>(sg[k*2], sg[k*2+1]);
						s += a*b;
						break;
					}
				}
			}

			sz = ijg[l+1] - ijg[l];

			if (sz == 1)
				a = std::complex<double>(gg[ijg[l]], 0);
			else
				a = std::complex<double>(gg[ijg[l]], gg[ijg[l]+1]);

			b = std::complex<double>(d[jg[l]*2], d[jg[l]*2+1]);

			temp = (a - s)/b;
			sg[l*2] = real(temp);
			sg[l*2+1] = imag(temp);
		}

		s = std::complex<double>(0, 0);
		for(k=ig[i]; k<=ig[i+1]-1; k++)
		{
			a = std::complex<double>(sg[k*2], sg[k*2+1]);
			s += a*a;
		}

		sz = idi[i+1] - idi[i];
		if (sz == 1)
			a = std::complex<double>(di[idi[i]], 0);
		else
			a = std::complex<double>(di[idi[i]], di[idi[i]+1]);

		temp = sqrt(a - s);

		d[i*2] = real(temp);
		d[i*2+1] = imag(temp);
	}

	cout << "Complex LLT end.\n";
	logfile << "Complex LLT end.\n";

	return 0;
}
//------------------------------------------------------------------------
// A forward solve with complex-valued sparse triangular matrix
//------------------------------------------------------------------------
int Block_2x2_solver::SolveL_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
				 double *f, double *x)
{
	int i, j, k;
	std::complex<double> s, _x, _t, _d, _f;

	for (i=0; i<nb; i++)
	{
		s = std::complex<double>(0, 0);

		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			_x = std::complex<double>(x[k*2], x[k*2+1]);
			_t = std::complex<double>(gg[j*2], gg[j*2+1]);
			s += _x * _t;
		}

		_f = std::complex<double>(f[i*2], f[i*2+1]);
		_d = std::complex<double>(di[i*2], di[i*2+1]);
		_x = (_f - s)/_d;
		x[i*2] = real(_x);
		x[i*2+1] = imag(_x);
	}	

	return 0;
}
//------------------------------------------------------------------------
// A backward solve with complex-valued sparse triangular matrix
//------------------------------------------------------------------------
int Block_2x2_solver::SolveU_Cmplx(int nb, int *ig, int *jg, double *di, double *gg,
				 double *f, double *s, double *x)
{
	int i, j, k;
	std::complex<double> _s, _x, _t, _d, _f;

	for(i=0; i<nb*2; i++)
		s[i] = 0;

	for(i=nb-1; i>=0; i--)
	{
		_f = std::complex<double>(f[i*2], f[i*2+1]);
		_s = std::complex<double>(s[i*2], s[i*2+1]);
		_d = std::complex<double>(di[i*2], di[i*2+1]);
		_x = (_f - _s)/_d;
		x[i*2] = real(_x);
		x[i*2+1] = imag(_x);

		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			_t = std::complex<double>(gg[j*2], gg[j*2+1]);
			_s = _x*_t;
			s[k*2] += real(_s);
			s[k*2+1] += imag(_s);
		}
	}

	return 0;
}
//------------------------------------------------------------------------
// Multiplication of upper triangular system with block-diagonal matrix by a vector
//------------------------------------------------------------------------
void Block_2x2_solver::mult_u_blockdiag(int nb, double *df, double *ggu_f, double *x, double *y)
{
	int i;
	int ind;

#pragma omp parallel shared(df, ggu_f, x, y)
	{
#pragma omp for private(i, ind)
		for(i=0; i<nb; i++)
		{
			ind = i*2;
			y[ind] = df[ind]*x[ind] + ggu_f[i]*x[ind+1];
			y[ind+1] = df[ind+1]*x[ind+1];
		}	
	} // end parallel section
}
//------------------------------------------------------------------------
// Multiplication of lower triangular system with block-diagonal matrix by a vector
//------------------------------------------------------------------------
void Block_2x2_solver::mult_l_blockdiag(int nb, double *df, double *ggl_f, double *x, double *y)
{
	int i;
	int ind;

#pragma omp parallel shared(df, ggl_f, x, y)
	{
#pragma omp for private(i, ind)
		for(i=0; i<nb; i++)
		{
			ind = i*2;
			y[ind]  =  df[ind]*x[ind];
			y[ind+1] = ggl_f[i]*x[ind] + df[ind+1]*x[ind+1];
		}	
	} // end parallel section
}
//------------------------------------------------------------------------
// Permutation of a vector
//------------------------------------------------------------------------
void Block_2x2_solver::Perm(int nb, double *x, double *y)
{
	for (int i=0; i<nb; i++)
	{
		y[i*2]   = -x[i*2+1];
		y[i*2+1] =  x[i*2];
	}
}
