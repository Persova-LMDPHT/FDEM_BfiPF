/*
 * GENERAL REMARKS
 *
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer. 
 *
 *                              FDEM_BfiPF
 *  The file contains code of functions of the class "Bound_cond"
 *  for imposing the boundary conditions in the 2D FEM harmonic global matrix and right hand side vector
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin                    
 *  Novosibirsk State Technical University,              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 * Version 2.0 December 10, 2024            
*/

#include "stdafx.h"
#include "bound_cond_2d.h"

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
//------------------------------------------------------ 
// Taking into account the first boundary conditions
//------------------------------------------------------ 
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

void Bound_Cond::Set_bound_cond_harm(int npls)
{
	int i, j, node, ipls;

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
			for(ipls=0;ipls<npls;ipls++)
			{
				pr[str*2+ipls*n_nodes_c*6] = 0.0;
				pr[str*2+1+ipls*n_nodes_c*6] = 0.0;
			}

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
//---------------------------------------------------------
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

void Bound_Cond::Set2kr(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
						double *mu2d, double current, int n_2kr, int (*edges_2kr)[2])
{
	int i;
	int v[2];
	double r;
	double f[2];
	double z[2];
	double coef;
	double hz;
	double val[2];

	for (i=0; i<n_2kr; i++)
	{
		v[0] = edges_2kr[i][0]-1;
		v[1] = edges_2kr[i][1]-1;

		r = xy[v[0]][0];

		z[0] = xy[v[0]][1];
		z[1] = xy[v[1]][1];

		hz = fabs(z[1] - z[0]);

		coef = (hz/6.0)/(2.0*PI);

		f[0] = 2.0*current*coef + 1.0*current*coef;
		f[1] = 1.0*current*coef + 2.0*current*coef;

		pr[v[0]] += f[0];
		pr[v[1]] += f[1];
	}
}

void Bound_Cond::Set2krHarm(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
							double *mu2d, double current, int n_2kr, int (*edges_2kr)[2],int ipls)
{
	int i;
	int v[2];
	double r;
	double f[2];
	double z[2];
	double coef;
	double hz;
	double val[2];

	for (i=0; i<n_2kr; i++)
	{
		v[0] = edges_2kr[i][0];
		v[1] = edges_2kr[i][1];

		r = xy[v[0]][0];

		z[0] = xy[v[0]][1];
		z[1] = xy[v[1]][1];

		hz = fabs(z[1] - z[0]);

		coef = (hz/6.0)/(2.0*PI);

		f[0] = 2.0*current*coef + 1.0*current*coef;
		f[1] = 1.0*current*coef + 2.0*current*coef;

		pr[v[0]*6+2+ipls*n_nodes_c*6] += f[0];
		pr[v[1]*6+2+ipls*n_nodes_c*6] += f[1];
	}
}

void Bound_Cond::Set2kr_3x3(int n_elem, int (*nvtr)[4], double (*xy)[2], int n_nodes_c,
						double *mu2d, double current, int n_2kr, int (*edges_2kr)[2])
{
	int i;
	int v[2];
	double r;
	double f[2];
	double z[2];
	double coef;
	double hz;
	double val[2];

	for (i=0; i<n_2kr; i++)
	{
		v[0] = edges_2kr[i][0]-1;
		v[1] = edges_2kr[i][1]-1;

		r = xy[v[0]][0];

		z[0] = xy[v[0]][1];
		z[1] = xy[v[1]][1];

		hz = fabs(z[1] - z[0]);

		coef = (hz/6.0)/(2.0*PI);

		f[0] = 2.0*current*coef + 1.0*current*coef;
		f[1] = 1.0*current*coef + 2.0*current*coef;

		pr[v[0]*3+2] += f[0];
		pr[v[1]*3+2] += f[1];
	}
}

void Bound_Cond::Set1kr(double (*xy)[2])
{

	int i, j, k;
	int node, str;
	double u_g;
	double *dirichlet=NULL;
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

	dirichlet = new double[n_nodes_c];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_bound_nodes; i++)
	{
		node = bound_nodes[i];
		dirichlet[node] = 0.0;
	}

	for(i=0; i<n_bound_nodes; i++)  
	{
		node = bound_nodes[i];
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

	for(i=0; i<n_nodes_c; i++)
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

void Bound_Cond::Set1kr_big(double (*xy)[2], double big)
{
	int i, j;
	double u_g;
	double *dirichlet=NULL;

	dirichlet = new double[n_nodes_c];
	if(dirichlet == 0)
		Memory_allocation_error("dirichlet", "Bound_Cond::Set_bound_cond_for_Hy_problem_with_symmetrization");

	for(i=0; i<n_bound_nodes; i++)
	{
		int node = bound_nodes[i];
		dirichlet[node] = 0.0;
	}

	for(i=0; i<n_bound_nodes; i++)  
	{
		int str = bound_nodes[i];
		u_g = dirichlet[str];

		di[str] = big;
		pr[str] = u_g*big;
	}

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

void Bound_Cond::SetDeltaFunction2(int node, double current, int (*nvtr)[4], double (*xy)[2], double omega)
{
	pr[node*6+4] -= current/(2.0*PI*omega);	
}

void Bound_Cond::SetDeltaFunctionElliptic(int node, double current, int (*nvtr)[4], double (*xy)[2], double omega)
{
	pr[node*6+2] += current/(2.0*PI);	
}

void Bound_Cond::SetDeltaFunction(int elem1, int elem2, double current, int (*nvtr)[4], double (*xy)[2], double omega)
{
	double r11, r12, r21, r22;
	double z11, z12, z21, z22;
	double hr1, hz1, hr2, hz2;
	double val;
	double coef1, coef2;

	r11 = xy[nvtr[elem1][0]][0];
	r12 = xy[nvtr[elem1][3]][0];
	z11 = xy[nvtr[elem1][0]][1];
	z12 = xy[nvtr[elem1][3]][1];

	r21 = xy[nvtr[elem2][0]][0];
	r22 = xy[nvtr[elem2][3]][0];
	z21 = xy[nvtr[elem2][0]][1];
	z22 = xy[nvtr[elem2][3]][1];

	hr1 = fabs(r12 - r11);
	hz1 = fabs(z12 - z11);

	hr2 = fabs(r22 - r21);
	hz2 = fabs(z22 - z21);

	val = current/(omega*PI*((r12*r12 - r11*r11)*hz1 + (r22*r22 - r21*r21)*hz2));

	SetF(elem1, val, nvtr, xy);
	SetF(elem2, val, nvtr, xy);
}

void Bound_Cond::SetF(int elem, double current, int (*nvtr)[4], double (*xy)[2])
{
	int i;
	double z1, z2, hz;
	double coef;

	z1 = xy[nvtr[elem][0]][1];
	z2 = xy[nvtr[elem][3]][1];
	hz = fabs(z2 - z1);

	cout << "hz=" << hz << endl;

	coef = (hz/6.0)/(2.0*PI);

	for (i=0; i<4; i++)
	{
		pr[nvtr[elem][i]*6+4] += 3.0*current*coef;
	}	
}
