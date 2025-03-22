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
 *  The file contains code of functions of the class "Bound_cond"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin          
 *  Novosibirsk State Technical University,                                                      
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                            
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                         
 *  Version 2.0 December 10, 2024                                                                 
*/

#include "stdafx.h"
#include "bound_cond_2d.h"
#include "in_out.h"
#include "For_Solvers.h"

double modul(complex<double> val)
{
	double re=real(val);
	double im=imag(val);
	return sqrt(re*re+im*im);
}

Bound_Cond::Bound_Cond(double *di, double *ggl, double *ggu, double *pr,
					   int *bound_nodes, int n_nodes_c, int n_bound_nodes,
					   int *ig, int *jg)
{
	this->di = di;
	this->ggl = ggl;
	this->ggu = ggu;
	this->pr = pr;
	this->bound_nodes = bound_nodes;
	this->n_nodes_c = n_nodes_c;
	this->n_bound_nodes = n_bound_nodes;
	this->ig = ig;
	this->jg = jg;
}

Bound_Cond::~Bound_Cond()
{
}

void Bound_Cond::Set_bound_cond_harm(double (*xy)[2])
{
	int i, j, node;
	double r;

	bool *is_node_bound=NULL;

	is_node_bound = new bool[n_nodes_c];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "Bound_Cond::Set_bound_cond_for_Ay_problem");

	for(i=0; i<n_nodes_c; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
	{
		node = bound_nodes[i];
		r = xy[node][0];
		if (r>0.001)
			is_node_bound[node] = true;
	}

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];
		r = xy[node][0];
		if (r<=0.001)
			continue;

		for(j=0; j<3; j++)
		{
			int str;
			str = node*3+j;

			di[str*2] = 1.0;
			di[str*2+1] = 0.0;
			pr[str*2] = 0.0;
			pr[str*2+1] = 0.0;

			for(int i1=ig[str]; i1<=ig[str+1]-1; i1++)
			{
				ggl[i1*2] = 0.0;
				ggl[i1*2+1] = 0.0;
			}
		}
	}

	for(int i1=0; i1<n_nodes_c*3; i1++)
	{
		for(int j1=ig[i1]; j1<=ig[i1+1]-1; j1++)
		{
			if(is_node_bound[jg[j1]/3]==true)
			{
				node = jg[j1]/3;
				if (r<=0.001)
					continue;

				ggu[j1*2] = 0.0;
				ggu[j1*2+1] = 0.0;
			}
		}
	}

	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
//------------------------------------------------------
// Taking into account the first boundary conditions
//------------------------------------------------------ 
void Bound_Cond::Set_bound_cond_harm_output(double (*u_bnd)[2],int npls)
{
	int i, j, k, ipls;
	int node, str;
	double u_re, u_im;
	double (*dirichlet)[2]=NULL;
	bool *is_dof_bound=NULL;

	is_dof_bound = new bool[n_nodes_c];
	if(is_dof_bound == 0)
		Memory_allocation_error("is_dof_bound", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c; i++)
		is_dof_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
	{
		is_dof_bound[bound_nodes[i]] = true;
	}

	dirichlet = new double[n_nodes_c][2];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_bound_nodes; i++)
	{
		node = bound_nodes[i];
		dirichlet[node][0] = u_bnd[i][0];
		dirichlet[node][1] = u_bnd[i][1];
	}

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];

		str = node;
		u_re = dirichlet[str][0];
		u_im = dirichlet[str][1];

		di[str*2]   = 1.0;
		di[str*2+1] = 0.0;

		for(ipls=0;ipls<npls;ipls++)
		{
			pr[str*2+n_nodes_c*ipls*2]   = u_re;
			pr[str*2+1+n_nodes_c*ipls*2] = u_im;
		}

		for(k=ig[str]; k<=ig[str+1]-1; k++)
		{
			if(is_dof_bound[jg[k]]==false)
			{
				for(ipls=0;ipls<npls;ipls++)
				{
					pr[jg[k]*2+n_nodes_c*ipls*2]   -= ggu[k*2]*u_re - ggu[k*2+1]*u_im;
					pr[jg[k]*2+1+n_nodes_c*ipls*2] -= ggu[k*2]*u_im + ggu[k*2+1]*u_re;
				}
				ggu[k*2]   = 0.0;
				ggu[k*2+1] = 0.0;
			}
			ggl[k*2]   = 0.0;
			ggl[k*2+1] = 0.0;
		}
	}

	for(i=0; i<n_nodes_c; i++)
	{
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if(is_dof_bound[k]==true) 
			{
				if(is_dof_bound[i] == false)
				{
					for(ipls=0;ipls<npls;ipls++)
					{
						pr[i*2+n_nodes_c*ipls*2]   -= ggl[j*2]*dirichlet[k][0] - ggl[j*2+1]*dirichlet[k][1];
						pr[i*2+1+n_nodes_c*ipls*2] -= ggl[j*2]*dirichlet[k][1] + ggl[j*2+1]*dirichlet[k][0];
					}
					ggl[j*2] = 0.0;
					ggl[j*2+1] = 0.0;
				}
				ggu[j*2] = 0.0;
				ggu[j*2+1] = 0.0;
			}
		}
	}

	if(is_dof_bound) {delete [] is_dof_bound; is_dof_bound=NULL;}
	if(dirichlet) {delete [] dirichlet; dirichlet=NULL;}
}

void Bound_Cond::Set_bound_cond_harm()
{
	int i, j, node;

	bool *is_node_bound=NULL;

	is_node_bound = new bool[n_nodes_c];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "Bound_Cond::Set_bound_cond_for_Ay_problem");

	for(i=0; i<n_nodes_c; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = true;

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];

		for(j=0; j<3; j++)
		{
			int str;
			str = node*3+j;

			di[str*2] = 1.0;
			di[str*2+1] = 0.0;
			pr[str*2] = 0.0;
			pr[str*2+1] = 0.0;

			for(int i1=ig[str]; i1<=ig[str+1]-1; i1++)
			{
				ggl[i1*2] = 0.0;
				ggl[i1*2+1] = 0.0;
			}
		}
	}

	for(int i1=0; i1<n_nodes_c*3; i1++)
	{
		for(int j1=ig[i1]; j1<=ig[i1+1]-1; j1++)
		{
			if(is_node_bound[jg[j1]/3]==true)
			{
				ggu[j1*2] = 0.0;
				ggu[j1*2+1] = 0.0;
			}
		}
	}

	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}
//------------------------------------------------------------------------
void Bound_Cond::Set_bound_cond_harm2(int n_2kr, int (*edges_2kr)[2])
{
	int i, j, node;
	bool *is_node_bound=NULL;
	int nBoundNodes;
	vector<int> boundNodes;

	is_node_bound = new bool[n_nodes_c];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "Bound_Cond::Set_bound_cond_for_Ay_problem");

	for(i=0; i<n_nodes_c; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_2kr; i++)
	{
		is_node_bound[edges_2kr[i][0]] = true;
		is_node_bound[edges_2kr[i][1]] = true;
	}

	nBoundNodes = 0;
	for(i=0; i<n_nodes_c; i++)
	{
		if (is_node_bound[i])
			nBoundNodes++;
	}
	boundNodes.resize(nBoundNodes);

	j = 0;
	for(i=0; i<n_nodes_c; i++)
	{
		if (is_node_bound[i])
		{
			boundNodes[j] = i;
			j++;
		}
	}


	for(i=0; i<nBoundNodes; i++)  
	{
		node = boundNodes[i];

		for(j=0; j<3; j++)
		{
			int str;
			str = node*3+j;

			if (j == 1)
				continue;

			di[str*2] = 1.0;
			di[str*2+1] = 0.0;
			pr[str*2] = 0.0;
			pr[str*2+1] = 0.0;

			for(int i1=ig[str]; i1<=ig[str+1]-1; i1++)
			{
				ggl[i1*2] = 0.0;
				ggl[i1*2+1] = 0.0;
			}
		}
	}

	for(int i1=0; i1<n_nodes_c*3; i1++)
	{
		for(int j1=ig[i1]; j1<=ig[i1+1]-1; j1++)
		{
			if(is_node_bound[jg[j1]/3]==true && jg[j1]%3 != 1)
			{
				ggu[j1*2] = 0.0;
				ggu[j1*2+1] = 0.0;
			}
		}
	}

	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}

void Bound_Cond::Set_bound_cond()
{
	int i, j, node;

	bool *is_node_bound=NULL;

	is_node_bound = new bool[n_nodes_c];
	if(is_node_bound == 0)
		Memory_allocation_error("is_node_bound", "Bound_Cond::Set_bound_cond_for_Ay_problem");

	for(i=0; i<n_nodes_c; i++)
		is_node_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
		is_node_bound[bound_nodes[i]] = true;

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];

		for(j=0; j<3; j++)
		{
			int str = node*3+j;
			di[str] = 1.0;
			pr[str] = 0.0;

			for(int i1=ig[str]; i1<=ig[str+1]-1; i1++)
			{
				ggl[i1] = 0.0;
			}
		}
	}

	for(int i1=0; i1<n_nodes_c*3; i1++)
	{
		for(int j1=ig[i1]; j1<=ig[i1+1]-1; j1++)
		{
			if(is_node_bound[jg[j1]/3]==true)
			{
				ggl[j1] = 0.0;
			}
		}
	}
	if(is_node_bound) {delete [] is_node_bound; is_node_bound=NULL;}
}

void Bound_Cond::Set_bound_cond_with_symmetrization(double (*xy)[2])
{
	int i, j, k;
	int node, str;
	double u_g;
	double *dirichlet=NULL;
	bool *is_dof_bound=NULL;

	is_dof_bound = new bool[n_nodes_c*3];
	if(is_dof_bound == 0)
		Memory_allocation_error("is_dof_bound", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c*3; i++)
		is_dof_bound[i] = false;

	for(i=0; i<n_bound_nodes; i++)
	{
		is_dof_bound[bound_nodes[i]*3] = true;
		is_dof_bound[bound_nodes[i]*3+1] = true;
		is_dof_bound[bound_nodes[i]*3+2] = true;
	}

	dirichlet = new double[n_nodes_c*3];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_bound_nodes; i++)
	{
		node = bound_nodes[i];
		dirichlet[node*3] = 0.0;
		dirichlet[node*3+1] = 0.0; 
		dirichlet[node*3+2] = 0.0; 
	}

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];
		for(j=0; j<3; j++)
		{
			str = node*3+j;
			u_g = dirichlet[str];

			di[str] = 1.0;
			pr[str] = u_g;

			for(k=ig[str]; k<=ig[str+1]-1; k++)
			{
				if(is_dof_bound[jg[k]]==false)
				{
					pr[jg[k]] -=ggu[k]*u_g;
					ggu[k] = 0.0;
				}
				ggl[k] = 0.0;
			}
		}
	}

	for(i=0; i<n_nodes_c*3; i++)
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if(is_dof_bound[k]==true) 
			{
				if(is_dof_bound[i] == false)
				{
					pr[i] -= ggl[j]*dirichlet[k];
					ggl[j] = 0.0;
				}
				ggu[j] = 0.0;
			}
		}	

	if(is_dof_bound) {delete [] is_dof_bound; is_dof_bound=NULL;}
	if(dirichlet) {delete [] dirichlet; dirichlet=NULL;}
}

void Bound_Cond::SwitchOffEquations(bool a1, bool a2, bool a3)
{
	int i, j, k;
	int node, str;
	double u_g;
	double *dirichlet=NULL;
	bool *is_dof_bound=NULL;

	is_dof_bound = new bool[n_nodes_c*3];
	if(is_dof_bound == 0)
		Memory_allocation_error("is_dof_bound", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c; i++)
	{
		is_dof_bound[i*3] = a1;
		is_dof_bound[i*3+1] = a2;
		is_dof_bound[i*3+2] = a3;
	}

	int n_bound_nodes_new = 0;
	if (a1)
		n_bound_nodes_new += n_nodes_c;
	if (a2)
		n_bound_nodes_new += n_nodes_c;
	if (a3)
		n_bound_nodes_new += n_nodes_c;
		
	dirichlet = new double[n_nodes_c*3];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c*3; i++)
		dirichlet[i] = 0.0;

	for(i=0; i<n_nodes_c*3; i++)  
	{
		if (!is_dof_bound[i])
			continue;

		node = i;
		str = node;
		u_g = dirichlet[str];

		di[str] = 1.0;
		pr[str] = u_g;

		for(k=ig[str]; k<=ig[str+1]-1; k++)
		{
			if(is_dof_bound[jg[k]]==false)
			{
				pr[jg[k]] -=ggu[k]*u_g;
				ggu[k] = 0.0;
			}
			ggl[k] = 0.0;
		}
	}

	for(i=0; i<n_nodes_c*3; i++)
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if(is_dof_bound[k]==true) 
			{
				if(is_dof_bound[i] == false)
				{
					pr[i] -= ggl[j]*dirichlet[k];
					ggl[j] = 0.0;
				}
				ggu[j] = 0.0;
			}
		}	

		if(is_dof_bound) {delete [] is_dof_bound; is_dof_bound=NULL;}
		if(dirichlet) {delete [] dirichlet; dirichlet=NULL;}
}

void Bound_Cond::SwitchOffEquationsHarm(bool a1, bool a2, bool a3)
{
	int i, j, k;
	int node, str;
	bool *is_dof_bound=NULL;

	is_dof_bound = new bool[n_nodes_c*3];
	if(is_dof_bound == 0)
		Memory_allocation_error("is_dof_bound", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_nodes_c; i++)
	{
		is_dof_bound[i*3] = a1;
		is_dof_bound[i*3+1] = a2;
		is_dof_bound[i*3+2] = a3;
	}

	int n_bound_nodes_new = 0;
	if (a1)
		n_bound_nodes_new += n_nodes_c;
	if (a2)
		n_bound_nodes_new += n_nodes_c;
	if (a3)
		n_bound_nodes_new += n_nodes_c;

	for(i=0; i<n_nodes_c*3; i++)  
	{
		if (!is_dof_bound[i])
			continue;

		node = i;
		str = node;

		di[str*2] = 1.0;
		di[str*2+1] = 0.0;
		pr[str*2] = 0.0;
		pr[str*2+1] = 0.0;

		for(k=ig[str]; k<=ig[str+1]-1; k++)
		{
			if(is_dof_bound[jg[k]]==false)
			{
				ggu[k*2] = 0.0;
				ggu[k*2+1] = 0.0;
			}
			ggl[k*2] = 0.0;
			ggl[k*2+1] = 0.0;
		}
	}

	for(i=0; i<n_nodes_c*3; i++)
		for(j=ig[i]; j<=ig[i+1]-1; j++)
		{
			k = jg[j];
			if(is_dof_bound[k]==true) 
			{
				if(is_dof_bound[i] == false)
				{
					ggl[j*2] = 0.0;
					ggl[j*2+1] = 0.0;
				}
				ggu[j*2] = 0.0;
				ggu[j*2+1] = 0.0;
			}
		}	

		if(is_dof_bound) {delete [] is_dof_bound; is_dof_bound=NULL;}
}

void Bound_Cond::CalcLocalMatrixUp(double (*m)[4], double hr, double hz, double rk, double gamma)
{
	double coef;
	double m11, m12, m22;

	coef = gamma*hr/(6.0*hz)/MU_0;

	m11 = coef * (2.0*rk + 0.5*hr);
	m12 = coef * (    rk + 0.5*hr);
	m22 = coef * (2.0*rk + 1.5*hr);

	m[0][0] = -m11;
	m[0][1] = -m12;
	m[0][2] =  m11;
	m[0][3] =  m12;

	m[1][0] = -m12;
	m[1][1] = -m22;
	m[1][2] =  m12;
	m[1][3] =  m22;

	m[2][0] = 0;
	m[2][1] = 0;
	m[2][2] = 0;
	m[2][3] = 0;
	
	m[3][0] = 0;
	m[3][1] = 0;
	m[3][2] = 0;
	m[3][3] = 0;
}

void Bound_Cond::CalcLocalMatrixDown(double (*m)[4], double hr, double hz, double rk, double gamma)
{
	double coef;
	double m11, m12, m22;

	coef = -gamma*hr/(6.0*hz)/MU_0;

	m11 = coef * (2.0*rk + 0.5*hr);
	m12 = coef * (    rk + 0.5*hr);
	m22 = coef * (2.0*rk + 1.5*hr);

	m[0][0] = 0;
	m[0][1] = 0;
	m[0][2] = 0;
	m[0][3] = 0;

	m[1][0] = 0;
	m[1][1] = 0;
	m[1][2] = 0;
	m[1][3] = 0;

	m[2][0] = m11;
	m[2][1] = m12;
	m[2][2] = -m11;
	m[2][3] = -m12;

	m[3][0] = m12;
	m[3][1] = m22;
	m[3][2] = -m12;
	m[3][3] = -m22;
}

void Bound_Cond::CalcLocalVectorUp(double *b, double hr, double rk, double gamma_re, double gamma_im)
{
	b[0] = gamma_re * 0.5*hr*(rk + hr/3.0)/MU_0;
	b[1] = gamma_im * 0.5*hr*(rk + hr/3.0)/MU_0;

	b[2] = gamma_re * 0.5*hr*(rk + hr*2.0/3.0)/MU_0;
	b[3] = gamma_im * 0.5*hr*(rk + hr*2.0/3.0)/MU_0;

	b[4] = 0;
	b[5] = 0;

	b[6] = 0;
	b[7] = 0;
}

void Bound_Cond::CalcLocalVectorDown(double *b, double hr, double rk, double gamma_re, double gamma_im)
{
	b[0] = 0;
	b[1] = 0;

	b[2] = 0;
	b[3] = 0;

	b[4] = gamma_re * 0.5*hr*(rk + hr/3.0)/MU_0;
	b[5] = gamma_im * 0.5*hr*(rk + hr/3.0)/MU_0;

	b[6] = gamma_re * 0.5*hr*(rk + hr*2.0/3.0)/MU_0;
	b[7] = gamma_im * 0.5*hr*(rk + hr*2.0/3.0)/MU_0;
}

void Bound_Cond::CalcLocalVectorErUp(double *b, double hr, double rk, double er_re, double er_im)
{
	b[0] = er_re * 0.5*hr*(rk + hr/3.0);
	b[1] = er_im * 0.5*hr*(rk + hr/3.0);

	b[2] = er_re * 0.5*hr*(rk + hr*(2.0/3.0));
	b[3] = er_im * 0.5*hr*(rk + hr*(2.0/3.0));

	b[4] = 0;
	b[5] = 0;

	b[6] = 0;
	b[7] = 0;
}

void Bound_Cond::CalcLocalVectorErDown(double *b, double hr, double rk, double er_re, double er_im)
{
	b[0] = 0;
	b[1] = 0;

	b[2] = 0;
	b[3] = 0;

	b[4] = er_re * 0.5*hr*(rk + hr/3.0);
	b[5] = er_im * 0.5*hr*(rk + hr/3.0);

	b[6] = er_re * 0.5*hr*(rk + hr*(2.0/3.0));
	b[7] = er_im * 0.5*hr*(rk + hr*(2.0/3.0));
}
//-----------------------------------------------------------
// Addition to the matrix and the right-hand side of the term with a jump in sigma at the separation of the media
//-----------------------------------------------------------
int Bound_Cond::SetSigma(int nelmJ, int *elmJ, double (*sigJ)[2], double (*epsJ)[2], int *lowUpElm, double k, int (*nvtr)[4], double (*xy)[2],
						 double *w_re, double *w_im, int npls,int kuzlov,int *nvkat,double *sigma2d,double *sigma2d_Z,double *dpr2d, double omega)
{
	int i, j, ipls;
	double m_up[4][4];
	double m_down[4][4];
	
	double *b_up=NULL;
	double *b_down=NULL;

	double gamma;
	int elm_up;
	int elm_down;
	complex<double> sig_up_r,sig_up_z;
	complex<double> sig_down_r,sig_down_z;

	double rk;
	double hr;
	double hz_up;
	double hz_down;
	complex<double> gamma_up;
	complex<double> gamma_down;

	double *dW_re=NULL;
	double *dW_im=NULL;

	int i_glob;
	int j_glob;

	int flag;
	int elm_cur;

	double ves;

	dW_re=new double[npls];
	dW_im=new double[npls];

	b_up=new double[8*npls];
	b_down=new double[8*npls];

	for(i=0;i<nelmJ/2;i++)
	{
		if(lowUpElm[i*2]!=0 || lowUpElm[i*2+1]!=1)
		{
			cout << "error"<<endl;
			cout << "nelmJ= "<<nelmJ<<endl;
			cout << "i= "<<i+1<<endl;
			exit(1);
		}

		elm_down = elmJ[i*2];
		elm_up   = elmJ[i*2+1];

		sig_down_r = complex<double> (sigma2d[nvkat[elm_down]],dpr2d[nvkat[elm_down]]*omega);
		sig_down_z = complex<double> (sigma2d_Z[nvkat[elm_down]],dpr2d[nvkat[elm_down]]*omega);
		sig_up_r   = complex<double> (sigma2d[nvkat[elm_up]],dpr2d[nvkat[elm_up]]*omega);
		sig_up_z   = complex<double> (sigma2d_Z[nvkat[elm_up]],dpr2d[nvkat[elm_up]]*omega);

		rk = xy[nvtr[elm_up][0]][0];
		hr = xy[nvtr[elm_up][1]][0] - rk;

		hz_up = xy[nvtr[elm_up][3]][1] - xy[nvtr[elm_up][0]][1];
		hz_down = xy[nvtr[elm_down][3]][1] - xy[nvtr[elm_down][0]][1];

		for(ipls=0;ipls<npls;ipls++)
		{
			dW_re[ipls] = (w_re[nvtr[elm_up][1]+ipls*kuzlov] - w_re[nvtr[elm_up][0]+ipls*kuzlov])/hr;
			dW_im[ipls] = (w_im[nvtr[elm_up][1]+ipls*kuzlov] - w_im[nvtr[elm_up][0]+ipls*kuzlov])/hr;
		}

		gamma_up = (sig_up_z - sig_down_z)/sig_up_r;
		gamma_down = (sig_up_z - sig_down_z)/sig_down_r;

		for (int i=0; i<4; i++)
		{
			for (j=0; j<4; j++)
			{
				m_up[i][j] = 0;
				m_down[i][j] = 0;
			}

			for(ipls=0;ipls<npls;ipls++)
			{
				b_down[i*2+ipls*8] = 0;
				b_down[i*2+1+ipls*8] = 0;

				b_up[i*2+ipls*8] = 0;
				b_up[i*2+1+ipls*8] = 0;
			}
		}

		if (modul(sig_up_r) < 1e-12 && modul(sig_down_r) > 1e-12)
		{
			CalcLocalMatrixDown(m_down, hr, hz_down, rk, 1.0);
			for(ipls=0;ipls<npls;ipls++)
			{
				complex<double> dW(dW_re[ipls],dW_im[ipls]);
				complex<double> MULT=dW*gamma_down;
				CalcLocalVectorDown(b_down+ipls*8, hr, rk, real(MULT), imag(MULT));
			}
			flag = -1;
		}
		else if (modul(sig_down_r) < 1e-12 && modul(sig_up_r) > 1e-12)
		{
			CalcLocalMatrixUp(m_up, hr, hz_up, rk, 1.0);
			for(ipls=0;ipls<npls;ipls++)
			{
				complex<double> dW(dW_re[ipls],dW_im[ipls]);
				complex<double> MULT=dW*gamma_up;
				CalcLocalVectorUp(b_up+ipls*8, hr, rk, real(MULT), imag(MULT));
			}
			flag = 1;	
		}
		else
		{
			CalcLocalMatrixDown(m_down, hr, hz_down, rk, 1.0);
			CalcLocalMatrixUp(m_up, hr, hz_up, rk, 1.0);

			for(ipls=0;ipls<npls;ipls++)
			{
				complex<double> dW(dW_re[ipls],dW_im[ipls]);
				complex<double> MULTDW=dW*gamma_down;
				complex<double> MULTUP=dW*gamma_up;

				CalcLocalVectorDown(b_down+ipls*8, hr, rk, real(MULTDW), imag(MULTDW));
				CalcLocalVectorUp(b_up+ipls*8, hr, rk, real(MULTUP), imag(MULTUP));   
			}

			ves = (modul(sig_down_r)) / (modul(sig_up_r) + modul(sig_down_r));

			for (int i=0; i<4; i++)
			{
				for (j=0; j<4; j++)
				{
					m_down[i][j] *= 1.0-ves;
					m_up[i][j] *= ves;
				}

				for(ipls=0;ipls<npls;ipls++)
				{
					b_down[i*2+ipls*8]   *= 1.0-ves;
					b_down[i*2+1+ipls*8] *= 1.0-ves;
					b_up[i*2+ipls*8]     *= ves;
					b_up[i*2+1+ipls*8]   *= ves;
				}
			}

			flag = 0;
		}

		switch (flag)
		{
		case -1:
			AddLocalMatrix(m_down, nvtr, elm_down, real(gamma_down), imag(gamma_down));
			for(ipls=0;ipls<npls;ipls++)
			{
				AddLocalVector(b_down+ipls*8, nvtr, elm_down, ipls);
			}
			break;
		case 1:
			AddLocalMatrix(m_up, nvtr, elm_up, real(gamma_up), imag(gamma_up));
			for(ipls=0;ipls<npls;ipls++)
			{
				AddLocalVector(b_up+ipls*8, nvtr, elm_up, ipls);
			}
			break;
		case 0:
			AddLocalMatrix(m_down, nvtr, elm_down, real(gamma_down), imag(gamma_down));
			AddLocalMatrix(m_up, nvtr, elm_up, real(gamma_up), imag(gamma_up));

			for(ipls=0;ipls<npls;ipls++)
			{
				AddLocalVector(b_down+ipls*8, nvtr, elm_down, ipls);
				AddLocalVector(b_up+ipls*8, nvtr, elm_up, ipls);
			}
			break;
		}
	}

	if(b_up){delete [] b_up; b_up=NULL;}
	if(b_down){delete [] b_down; b_down=NULL;}

	if(dW_re){delete [] dW_re; dW_re=NULL;}
	if(dW_im){delete [] dW_im; dW_im=NULL;}

	return 0;
}

void Bound_Cond::AddLocalMatrix(double (*m)[4], int (*nvtr)[4], int elm, double gamma_re,double gamma_im)
{
	int i, j;
	int i_glob, j_glob;

	for (i=0; i<4; i++)
	{
		i_glob = nvtr[elm][i];

		di[i_glob*2] += m[i][i]*gamma_re; 
		di[i_glob*2+1] += m[i][i]*gamma_im; 

		for (j=0; j<4; j++)
		{
			j_glob = nvtr[elm][j];

			if(j_glob < i_glob)
			{
				for(int mm=ig[i_glob]; mm<=ig[i_glob+1]-1; mm++)
				{
					if(jg[mm]==j_glob)
					{
						ggl[mm*2] += m[i][j]*gamma_re;
						ggl[mm*2+1] += m[i][j]*gamma_im;
						ggu[mm*2] += m[j][i]*gamma_re;
						ggu[mm*2+1] += m[j][i]*gamma_im;
					}
				}
			}
		}
	}
}

void Bound_Cond::AddLocalVector(double *v, int (*nvtr)[4], int elm, int ipls)
{
	int i;
	int i_glob;

	for (i=0; i<4; i++)
	{
		i_glob = nvtr[elm][i];

		pr[i_glob*2+ipls*2*n_nodes_c]   += v[i*2]; 
		pr[i_glob*2+1+ipls*2*n_nodes_c] += v[i*2+1]; 
	}
}
