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
 *  This file contains subroutines for the folded preconditioner       
 *                                                                                                                        
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/FoldedPreconditioner.cpp                                                                                     
 *  Novosibirsk State Technical University,                                                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                                    
 *  Version 1.2 April 9, 2021                                                                                             
*/                                                                                                                        

#include "stdafx.h"
#include "FoldedPreconditioner.h"
#include "ControlOMP.h"
extern ControlOMP omp;
extern ofstream logfile;
//------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------
FoldedPreconditioner::FoldedPreconditioner()
{
	n = 0;
	m = 0;

	zero = 1e-30;

	iptr_g = NULL;
	jptr_g = NULL;
	gg_g = NULL;
	iptr_gt = NULL;
	jptr_gt = NULL;
	gg_gt = NULL;
	diag_a = NULL;
	diag_gagt = NULL;
	help = NULL;
	help2 = NULL;
	help3 = NULL;
}
//------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------
FoldedPreconditioner::~FoldedPreconditioner()
{
	if (iptr_g) {delete [] iptr_g; iptr_g=NULL;}
	if (jptr_g) {delete [] jptr_g; jptr_g=NULL;}
	if (gg_g) {delete [] gg_g; gg_g=NULL;}
	if (iptr_gt) {delete [] iptr_gt; iptr_gt=NULL;}
	if (jptr_gt) {delete [] jptr_gt; jptr_gt=NULL;}
	if (gg_gt) {delete [] gg_gt; gg_gt=NULL;}
	if (diag_a) {delete [] diag_a; diag_a=NULL;}
	if (diag_gagt) {delete [] diag_gagt; diag_gagt=NULL;}
	if (help) {delete [] help; help=NULL;}
	if (help2) {delete [] help2; help2=NULL;}
	if (help3) {delete [] help3; help3=NULL;}
}
//------------------------------------------------------------------------
// Construction of the discrete gradient matrix
//------------------------------------------------------------------------
void FoldedPreconditioner::BuildGMatrix(int kpar, int n_edges_c, int n_nodes_c,
				   int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
				   int *ig_t, int *jg_t, double *gg_t, int *is_node_bound)
{
	int i, j, k;
	int v;
	int v_renum;
	vector<int> isNodeInSigma;

	if (iptr_g) {delete [] iptr_g; iptr_g=NULL;}
	if (jptr_g) {delete [] jptr_g; jptr_g=NULL;}
	if (gg_g) {delete [] gg_g; gg_g=NULL;}
	if (iptr_gt) {delete [] iptr_gt; iptr_gt=NULL;}
	if (jptr_gt) {delete [] jptr_gt; jptr_gt=NULL;}
	if (gg_gt) {delete [] gg_gt; gg_gt=NULL;}
	if (help) {delete [] help; help=NULL;}

	n = n_edges_c;

	// finding the number of (conforming) nodes in a conductive medium

	isNodeInSigma.resize(n_nodes_c, 0);

	for (i=0; i<kpar; i++)
	{
		int node;

		if (sigma3d[nvkat[i]] > 0.0)
		{
			for (j=0; j<8; j++)
			{
				node = nver[i][j];

				if (node < n_nodes_c)
				{
					// skip the nodes lying on the boundary of the computational domain
					if (is_node_bound[node])
						continue;

					isNodeInSigma[node] = 1;
				}
			}
		}
	}

	// the nodes lying in the conductive medium are numbered first
	vector<int> nodeRenum;
	nodeRenum.resize(n_nodes_c, -1);

	k=0;
	for (i=0; i<n_nodes_c; i++)
	{
		if (isNodeInSigma[i])
		{
			nodeRenum[i] = k;
			k++;
		}			
	}

	m = k; // save the number of nodes in the conductive medium

	// fill iptr_gt[]
	iptr_gt = new int[n_edges_c+1];
	iptr_gt[0] = 0;

	for (i=0; i<n_edges_c; i++)
	{
		j = 0;

		for (k=0; k<2; k++)
		{
			v = edges[i][k];

			if (v < n_nodes_c)
			{
				if(isNodeInSigma[v])
					j++;
			}
			else
			{
				for (int ii=ig_t[v]; ii<ig_t[v+1]; ii++)
				{
					if (isNodeInSigma[jg_t[ii]])
						j++;
				}
			}
		}

		iptr_gt[i+1] = iptr_gt[i] + j;
	}

	jptr_gt = new int[iptr_gt[n_edges_c]];
	gg_gt = new double[iptr_gt[n_edges_c]];

	// fill jptr_gt[] and gg_gt[]
	j = 0;
	for (i=0; i<n_edges_c; i++)
	{
		for (k=0; k<2; k++)
		{
			v = edges[i][k];

			if (v < n_nodes_c)
			{
				if (isNodeInSigma[v])
				{
					v_renum = nodeRenum[v];

					if (v_renum == -1)
						logfile << "FATAL ERROR IN G-MATRIX!!!!!!!\n" << flush;

					jptr_gt[j] = v_renum;

					if (k==0)
					{
						gg_gt[j] = -0.5;
					}
					else
					{	
						gg_gt[j] = 0.5;
					}
					j++;
				}
			}
			else
			{
				for (int ii=ig_t[v]; ii<ig_t[v+1]; ii++)
				{
					int kk = jg_t[ii];

					if (isNodeInSigma[kk])
					{
						v_renum = nodeRenum[kk];

						if (v_renum == -1)
							logfile << "FATAL ERROR IN G-MATRIX!!!!!!!\n" << flush;

						jptr_gt[j] = v_renum;

						if (k == 0)
						{
							gg_gt[j] = -0.5*gg_t[ii];
						}
						else
						{
							gg_gt[j] = 0.5*gg_t[ii];
						}
						
						j++;
					}
				}
			}
		}
	}

	// count how many elements are in each column
	vector<int> temp;
	temp.resize(m, 0);

	for (i=0; i<n_edges_c; i++)
	{
		for (j=iptr_gt[i]; j<iptr_gt[i+1]; j++)
		{
			temp[jptr_gt[j]]++;
		}
	}

	// transpose
	iptr_g = new int[m+1];
	jptr_g = new int[iptr_gt[n_edges_c]];
	gg_g = new double[iptr_gt[n_edges_c]];

	iptr_g[0] = 0;
	for (i=0; i<m; i++)
		iptr_g[i+1] = iptr_g[i] + temp[i];

	for (i=0; i<m; i++)
		temp[i] = iptr_g[i];

	for (i=0; i<n_edges_c; i++)
	{
		for (j=iptr_gt[i]; j<iptr_gt[i+1]; j++)
		{
			k = jptr_gt[j];

			gg_g[temp[k]] = gg_gt[j];
			jptr_g[temp[k]] = i;

			temp[k]++;
		}
	}	
}
//------------------------------------------------------------------------
// / Build the main diagonal of the matrix P(G+M)PT
//------------------------------------------------------------------------
void FoldedPreconditioner::BuildDiag(int *ig, int *jg, double *di, double *gg)
{
	int i, j, k;
	int i2, j2;
	int i3;
	double sum;

	int *iptr_a = NULL;
	int *jptr_a = NULL;
	double *aelem = NULL;
	int sz_iptr;
	int sz_jptr;

	FromRSFtoCSR_1(n, ig, &sz_iptr, &sz_jptr);
	iptr_a = new int[sz_iptr];
	jptr_a = new int[sz_jptr];
	aelem = new double[sz_jptr];
	FromRSFtoCSR_2(n, ig, jg, di, gg, gg, iptr_a, jptr_a, aelem);



	if (help) {delete [] help; help=NULL;}
	if (diag_a) {delete [] diag_a; diag_a=NULL;}
	if (diag_gagt) {delete [] diag_gagt; diag_gagt=NULL;}

	help = new double[m];
	diag_a = new double[n];
	diag_gagt = new double[m];


	for (i=0; i<n; i++)
	{
		diag_a[i] = 1.0/di[i];
	}

	for (i=0; i<m; i++)
	{
		diag_gagt[i] = 0;
	}

	for (i=0; i<m; i++)
	{
		for (j=iptr_g[i]; j<iptr_g[i+1]; j++)
		{
			k = jptr_g[j];

			sum = 0;

			// multiply k-th row of matrix A by i-th column of matrix G^T
			for (i2=iptr_a[k]; i2<iptr_a[k+1]; i2++)
			{
				j2 = jptr_a[i2];

				// find element G(i, j2)
				for (i3=iptr_g[i]; i3<iptr_g[i+1]; i3++)
				{
					if (jptr_g[i3]==j2)
					{
						sum += aelem[i2]*gg_g[i3];
						break;
					}
				}
			}

			diag_gagt[i] += sum*gg_g[j];
		}
	}

	for (i=0; i<m; i++)
	{
		if (fabs(diag_gagt[i]) > zero)
		{
			diag_gagt[i] = 1.0/diag_gagt[i];
		}
		else
		{
			logfile << "diag_gagt[" << i << "] = " << diag_gagt[i] << endl;
			diag_gagt[i] = 1.0;
		}
	}

	if (iptr_a) {delete [] iptr_a; iptr_a=NULL;}
	if (jptr_a) {delete [] jptr_a; jptr_a=NULL;}
	if (aelem)  {delete [] aelem; aelem=NULL;}
}
//------------------------------------------------------------------------
// Applying the folded preconditioner to a vector
//------------------------------------------------------------------------
void FoldedPreconditioner::ApplyPreconditioner(double *x, double *y)
{
	int i, j, k;

#pragma omp parallel shared(x, y) private(j, k) num_threads(omp.GetNumberOfThreads())
{
	#pragma omp for nowait
	for (i=0; i<n; i++)
	{
		y[i] = x[i]*diag_a[i];
	}

	#pragma omp for
	for (i=0; i<m; i++)
	{
		help[i] = 0;
	}

	#pragma omp for
	for (i=0; i<m; i++)
	{
		for (j=iptr_g[i]; j<iptr_g[i+1]; j++)
		{
			k = jptr_g[j];
			help[i] += gg_g[j]*x[k];
		}
	}

	#pragma omp for
	for (i=0; i<m; i++)
	{
		help[i] *= diag_gagt[i];
	}

	#pragma omp for
	for (i=0; i<n; i++)
	{
		for (j=iptr_gt[i]; j<iptr_gt[i+1]; j++)
		{
			k = jptr_gt[j];
			y[i] += gg_gt[j]*help[k];
		}
	}
}// end parallel
}
//------------------------------------------------------------------------
// Conversion to the CSR format
//------------------------------------------------------------------------
void FoldedPreconditioner::FromRSFtoCSR_1(int n, int *ig, int *sz_iptr, int *sz_jptr)
{
	*sz_iptr = n+1;
	*sz_jptr = ig[n]*2 + n;
}
//------------------------------------------------------------------------
// Conversion to the CSR format
//------------------------------------------------------------------------
void FoldedPreconditioner::FromRSFtoCSR_2(int n, int *ig, int *jg, double *di, double *ggl, double *ggu, int *iptr, int *jptr, double *aelem)
{
	int i, j, k;
	vector<int> col; 

	// counting the number of elements in each row
	col.resize(n, 0);

	for (i=0; i<n; i++)
	{
		col[i] += ig[i+1] - ig[i] + 1; // lower triangle + diagonal

		// upper triangle
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			col[k]++;
		}
	}

	// iptr
	iptr[0] = 0;
	for (i=0; i<n; i++)
		iptr[i+1] = iptr[i] + col[i];

	if (iptr[n] != ig[n]*2 + n)
	{
		cout << "internal error: iptr[n]=" << iptr[n] << "ig[n]*2+n=" << ig[n]*2+n << endl;
	}

	// jptr, aelem
	for (i=0; i<n; i++)
		col[i] = iptr[i]; // in which position to put the value

	for (i=0; i<n; i++)
	{
		// lower triangle
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[i]] = k;
			aelem[col[i]] = ggl[j];
			col[i]++;
		}

		// diagonal
		jptr[col[i]] = i;
		aelem[col[i]] = di[i];
		col[i]++;
	}

	// upper triangle
	for (i=0; i<n; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[k]] = i;
			aelem[col[k]] = ggu[j];
			col[k]++;
		}
	}
}
//------------------------------------------------------------------------
// Build the main complex-valued diagonal of the matrix P(G+iwM)PT
//------------------------------------------------------------------------
void FoldedPreconditioner::BuildDiagComplex(int *ig, int *jg, int *idi, int *ijg, double *di, double *gg)
{
	int i, j, k;
	int i2, j2;
	int i3;
	std::complex<double> sum;
	std::complex<double> x, y, z;

	int *iptr_a = NULL;
	int *jptr_a = NULL;
	int *ijptr_a = NULL;
	double *aelem = NULL;
	int sz_iptr;
	int sz_jptr;
	int sz_ijptr;
	int sz_aelem;

	From2x2ToCSR2x2_1(n, ig, idi, ijg, &sz_iptr, &sz_jptr, &sz_ijptr, &sz_aelem);
	iptr_a = new int[sz_iptr];
	jptr_a = new int[sz_jptr];
	ijptr_a = new int[sz_ijptr];
	aelem = new double[sz_aelem];
	From2x2ToCSR2x2_2(n, ig, jg, idi, ijg, di, gg, iptr_a, jptr_a, ijptr_a, aelem);

	if (help) {delete [] help; help=NULL;}
	if (diag_a) {delete [] diag_a; diag_a=NULL;}
	if (diag_gagt) {delete [] diag_gagt; diag_gagt=NULL;}

	help = new double[2*m];
	diag_a = new double[2*n];
	diag_gagt = new double[2*m];



	for (i=0; i<n; i++)
	{
		int sz = idi[i+1] - idi[i];

		if (sz == 1)			
			x = std::complex<double>(di[idi[i]], 0);
		if (sz == 2)
			x = std::complex<double>(di[idi[i]], di[idi[i]+1]);

		y = 1.0/x;

		diag_a[i*2] = y.real();
		diag_a[i*2+1] = y.imag();
	}

	for (i=0; i<m*2; i++)
	{
		diag_gagt[i] = 0;
	}

	for (i=0; i<m; i++)
	{
		for (j=iptr_g[i]; j<iptr_g[i+1]; j++)
		{
			k = jptr_g[j];

			sum = 0;

			// multiply k-th row of matrix A by i-th column of matrix G^T
			for (i2=iptr_a[k]; i2<iptr_a[k+1]; i2++)
			{
				j2 = jptr_a[i2];

				// find element G(i, j2)
				for (i3=iptr_g[i]; i3<iptr_g[i+1]; i3++)
				{
					if (jptr_g[i3]==j2)
					{
						int sz = ijptr_a[i2+1] - ijptr_a[i2];

						if (sz == 1)
							x = std::complex<double>(aelem[ijptr_a[i2]], 0);
						if (sz == 2)
							x = std::complex<double>(aelem[ijptr_a[i2]], aelem[ijptr_a[i2]+1]);
						
						sum += x*gg_g[i3];
					}
				}
			}

			diag_gagt[i*2] += sum.real()*gg_g[j];
			diag_gagt[i*2+1] += sum.imag()*gg_g[j];
		}
	}

	for (i=0; i<m; i++)
	{
		x = std::complex<double>(diag_gagt[i*2], diag_gagt[i*2+1]);
		y = 1.0/x;
		diag_gagt[i*2] = y.real();
		diag_gagt[i*2+1] = y.imag();
	}

	if (iptr_a) {delete [] iptr_a; iptr_a=NULL;}
	if (jptr_a) {delete [] jptr_a; jptr_a=NULL;}
	if (ijptr_a) {delete [] ijptr_a; ijptr_a=NULL;}
	if (aelem)  {delete [] aelem; aelem=NULL;}
}
//------------------------------------------------------------------------
// Applying the complex-valued folded preconditioner to a vector
//------------------------------------------------------------------------
void FoldedPreconditioner::ApplyPreconditionerComplex(double *x, double *y)
{
	int i, j, k;
	double t1, t2;

#pragma omp parallel shared(x, y) private(j, t1, t2, k) num_threads(omp.GetNumberOfThreads())
{
	#pragma omp for nowait
	for (i=0; i<n; i++)
	{
		y[i*2] = x[i*2]*diag_a[i*2] - x[i*2+1]*diag_a[i*2+1];
		y[i*2+1] = x[i*2]*diag_a[i*2+1] + x[i*2+1]*diag_a[i*2];
	}

	#pragma omp for
	for (i=0; i<m*2; i++)
	{
		help[i] = 0;
	}

	#pragma omp for
	for (i=0; i<m; i++)
	{
		for (j=iptr_g[i]; j<iptr_g[i+1]; j++)
		{
			k = jptr_g[j];
			help[i*2]   += gg_g[j]*x[k*2];
			help[i*2+1] += gg_g[j]*x[k*2+1];
		}
	}

	#pragma omp for
	for (i=0; i<m; i++)
	{
		t1 = help[i*2]*diag_gagt[i*2] - help[i*2+1]*diag_gagt[i*2+1];
		t2 = help[i*2+1]*diag_gagt[i*2] + help[i*2]*diag_gagt[i*2+1];
		help[i*2] = t1;
		help[i*2+1] = t2;
	}

	#pragma omp for
	for (i=0; i<n; i++)
	{
		for (j=iptr_gt[i]; j<iptr_gt[i+1]; j++)
		{
			k = jptr_gt[j];
			y[i*2]   += gg_gt[j]*help[k*2];
			y[i*2+1] += gg_gt[j]*help[k*2+1];
		}
	}
}
}
//------------------------------------------------------------------------
// Conversion to the block CSR format
//------------------------------------------------------------------------
void FoldedPreconditioner::From2x2ToCSR2x2_1(int n, int *ig, int *idi, int *ijg,
										int *sz_iptr, int *sz_jptr, int *sz_ijptr, int *sz_aelem)
{
	*sz_iptr = n+1;
	*sz_jptr = ig[n]*2 + n;
	*sz_ijptr = ig[n]*2 + n + 1;
	*sz_aelem = ijg[ig[n]]*2 + idi[n];
}
//------------------------------------------------------------------------
// Conversion to the block CSR format
//------------------------------------------------------------------------
void FoldedPreconditioner::From2x2ToCSR2x2_2(int n, int *ig, int *jg, int *idi, int *ijg,
										double *di_block, double *gg_block, 
										int *iptr, int *jptr, int *ijptr, double *aelem)
{
	int sz;
	int i, j, k, m;
	vector<int> col; 
	vector<int> col2; 

	// counting the number of elements in each row
	col.resize(n, 0);
	col2.resize(n, 0); // rows sizes

	for (i=0; i<n; i++)
	{
		col[i] += ig[i+1] - ig[i] + 1; // lower triangle + diagonal
		col2[i] += (ijg[ig[i+1]] - ijg[ig[i]]) + (idi[i+1] - idi[i]);

		// upper triangle
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			col[k]++;
			col2[k] += ijg[j+1] - ijg[j];
		}
	}


	// iptr
	iptr[0] = 0;
	for (i=0; i<n; i++)
		iptr[i+1] = iptr[i] + col[i];

	if (iptr[n] != ig[n]*2 + n)
	{
		cout << "internal error: iptr[n]=" << iptr[n] << "ig[n]*2+n=" << ig[n]*2+n << endl;
	}


	// jptr, ijptr, aelem
	ijptr[0] = 0;
	for (i=0; i<n; i++)
		ijptr[iptr[i+1]] = ijptr[iptr[i]] + col2[i]; //in which position to put the value of aelem

	for (i=0; i<n; i++)
		col[i] = iptr[i]; // in which position to put the value

	ijptr[0] = 0;
	for (i=0; i<n; i++)
	{
		// lower triangle
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			jptr[col[i]] = jg[j];

			sz = ijg[j+1] - ijg[j];

			ijptr[col[i]+1] = ijptr[col[i]] + sz;

			for (k=0; k<sz; k++)
				aelem[ijptr[col[i]] + k] = gg_block[ijg[j] + k];

			col[i]++;
		}

		// diagonal
		jptr[col[i]] = i;
		sz = idi[i+1] - idi[i];
		ijptr[col[i]+1] = ijptr[col[i]] + sz;

		for (k=0; k<sz; k++)
			aelem[ijptr[col[i]] + k] = di_block[idi[i] + k];

		col[i]++;
	}

	// upper triangle
	for (i=0; i<n; i++)
	{
		for (j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			jptr[col[k]] = i;
			sz = ijg[j+1] - ijg[j];
			ijptr[col[k]+1] = ijptr[col[k]] + sz;

			for (m=0; m<sz; m++)
				aelem[ijptr[col[k]] + m] = gg_block[ijg[j] + m];

			col[k]++;
		}
	}	
}
//------------------------------------------------------------------------
// Initialization of the folded preconditioner
//------------------------------------------------------------------------
void FoldedPreconditioner::Prepare(int kpar, int n_edges_c, int n_nodes_c,
								   int *nvkat, int (*nver)[14], double *sigma3d, int (*edges)[2],
								   int *ig_t, int *jg_t, double *gg_t, int *is_node_bound,
								   int *ig, int *jg, double *di, double *gg)
{
	BuildGMatrix(kpar, n_edges_c, n_nodes_c, nvkat, nver, sigma3d, edges, ig_t, jg_t, gg_t, is_node_bound);
}
