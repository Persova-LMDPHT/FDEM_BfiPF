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
 *  This file contains code of functions for building subdomains for smoothing 3D EM field
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/Subdomain.cpp
 *  Modified by D.Sc. Denis V. Vagin 
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024
*/

#include "stdafx.h"
#include "Subdomain.h"
#include "Hex_Local_Matrix.h"
#include "PointVector.h"
#include "gauss_3.h"
#define EpsComp 1e-6

extern ofstream logfile;

int Subdomain::subind=0;
int Subdomain::npls=0;
int Subdomain::ntimes=0;

struct double_eps
{
	double value;
	void operator = (double &a){value=a;}
	double get_value(){return value;}
};

bool operator == (double_eps &a,double_eps &b){return fabs(a.value-b.value)<EpsComp;}
bool operator > (double_eps &a,double_eps &b){return a.value>b.value+EpsComp;}
bool operator < (double_eps &a,double_eps &b){return a.value<b.value-EpsComp;}
bool operator >= (double_eps &a,double_eps &b){return (a>b || a==b);}
bool operator <= (double_eps &a,double_eps &b){return (a<b || a==b);}

double_eps operator + (double_eps &a,double_eps &b)
{
	double_eps t;
	t.value=a.value+b.value;
	return t;
}

double_eps operator - (double_eps &a,double_eps &b)
{
	double_eps t;
	t.value=a.value-b.value;
	return t;
}

double operator * (double_eps &a,double &b)
{
	return a.value*b;
}

double operator * (double_eps &a,double_eps &b)
{
	return a.value*b.value;
}

bool SIGMANCOMP(_SIGMA_N el_1,_SIGMA_N el_2){return el_1.gtn < el_2.gtn;}
int Sn[6][4]={{0,2,4,6},{0,1,4,5},{0,1,2,3},{1,3,5,7},{2,3,6,7},{4,5,6,7}};

Subdomain::Subdomain()
{
	nver = NULL;
	xyz = NULL;
	xyz_r = NULL;
	p = NULL;
	di = NULL;
	gg = NULL;
	pr = NULL;
	x = NULL;
	d = NULL;
	sg = NULL;

	subind++;
	my_ind=subind;
}

Subdomain::~Subdomain()
{
	if (p)    {delete p; p = NULL;}
	if (pr)   {delete [] pr; pr = NULL;}
	if (x)   {delete [] x; x = NULL;}
	if (di)   {delete [] di; di = NULL;}
	if (gg)   {delete [] gg; gg = NULL;}
	
	if (xyz)  {delete [] xyz;  xyz  = NULL;}
	if (xyz_r){delete [] xyz_r;xyz_r= NULL;}
	if (nver) {delete [] nver; nver = NULL;}

	if (d)    {delete [] d; d = NULL;}
	if (sg)   {delete [] sg; sg = NULL;}
}

int Subdomain::Init(int material, int levelNeighbors, vector< vector<long> > &PointresForElem, AbstractFEM3D *TaskCalcMesh)
{
	int i, j, k, m, l, rr;
	int level;
	bool flag;
	int sz;
	vector<long> newElems, newElems2;
	vector<bool> isElemInSubdomain;
	vector<bool> isNodeInSubdomain;

	this->material = material;

	int emat;
	int nnoe=TaskCalcMesh->GetElementNodesNumber();
	int kpar=TaskCalcMesh->GetNumberOfElements();
	int kuzlov=TaskCalcMesh->GetNumberOfNodes();

	isElemInSubdomain.resize(kpar, false);
	isNodeInSubdomain.resize(kuzlov, false);

	vector<int> renumNodeFromOldToNew;
	renumNodeFromOldToNew.resize(kuzlov, -1);

	for (i=0; i<kpar; i++)
	{
		if (PointresForElem[i].size() > 0)
		{
			emat=TaskCalcMesh->GetElementMaterial(i);

			if (this->material != -1)
			{
				if (this->material != emat)
					continue;
			}

			isElemInSubdomain[i] = true;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(i,j);
				if (k >= 0) 
					isNodeInSubdomain[k] = true;
			}
		}
	}

	for (level=0; level < levelNeighbors; level++)
	{
		newElems.clear();
		newElems2.clear();

		for (i=0; i<kpar; i++)
		{
			if (isElemInSubdomain[i])
				continue;

			if (this->material != -1)
			{
				emat=TaskCalcMesh->GetElementMaterial(i);

				if (this->material != emat)
					continue;
			}

			flag = false;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(i,j);

				if (k >= 0)
				{
					if (isNodeInSubdomain[k])
					{
						flag = true;
						break;
					}
				}
			}

			if (flag)
				newElems.push_back(i);
		}

		std::sort(newElems.begin(), newElems.end());
		std::unique_copy(newElems.begin(), newElems.end(), back_inserter(newElems2));

		sz = (int)newElems2.size();
		for (i=0; i<sz; i++)
		{
			isElemInSubdomain[newElems2[i]] = true;

			for (j=0; j<nnoe; j++)
			{
				k = TaskCalcMesh->GetNodeNumberOnElement(newElems2[i],j);

				if (k >= 0)
					isNodeInSubdomain[k] = true;
			}
		}
	}

	this->n_elem = 0;
	for (i=0; i<kpar; i++)
	{
		if (isElemInSubdomain[i])
			this->n_elem++;
	}

	this->ValueInCenter.resize(this->n_elem);
	this->renumElemFromNewToOld.resize(this->n_elem, -1);

	this->n_nodes = 0;
	for (i=0; i<kuzlov; i++)
	{
		if (isNodeInSubdomain[i])
			this->n_nodes++;
	}

	this->renumNodeFromNewToOld.resize(this->n_nodes, -1);

	k = 0;
	for (i=0; i<kpar; i++)
	{
		if (isElemInSubdomain[i])
		{
			this->renumElemFromNewToOld[k] = i;
			k++;
		}
	}

	renumElemFromOldToNew.clear();
	renumElemFromOldToNew.resize(kpar,-1);

	for (i=0; i<this->n_elem; i++)renumElemFromOldToNew[renumElemFromNewToOld[i]]=i;
	
	if (this->xyz)  {delete [] this->xyz;  this->xyz = NULL;}
	this->xyz = new double[kuzlov][3];

	if (this->xyz_r)  {delete [] this->xyz_r;  this->xyz_r = NULL;}
	this->xyz_r = new double[kuzlov][3];

	if (this->nver) {delete [] this->nver; this->nver = NULL;}
	this->nver = new long[kpar][14];

	vector<_SIGMA_N> PRE_SIGMA_N,SIGMA_N;

	PRE_SIGMA_N.clear();
	SIGMA_N.clear();

	ElemNeibVec.clear();
	ElemNeibVec.resize(n_elem);


	vector<ElemNeib> ElemNeibVecTmp;
	ReadElemNeib(ElemNeibVecTmp,kpar);
	for(i=0;i<kpar;i++){
		if(isElemInSubdomain[i]){
			m=renumElemFromOldToNew[i];
			for(j=0;j<6;j++){
				sz=(int)ElemNeibVecTmp[i].neib[j].size();
				for(k=0;k<sz;k++){
					l=ElemNeibVecTmp[i].neib[j][k];
					if(l<0 || l>=(int)isElemInSubdomain.size())
					{
						cout<<i+1<<"---"<<l+1<<endl;
						cout<<(int)isElemInSubdomain.size()<<endl;
						cout<<"stop"<<endl;
						exit(1001);
					}

					if(isElemInSubdomain[l]){
						ElemNeibVec[m].neib[j].push_back(renumElemFromOldToNew[l]);
					}
				}
			}
		}
	}

	for (i=0; i<kpar; i++){
		for (j=0; j<nnoe; j++){
			pv::Point3D TempPoint;
			k = TaskCalcMesh->GetNodeNumberOnElement(i,j);
			nver[i][j]=k;
			TempPoint=TaskCalcMesh->GetNode(k);
			this->xyz[k][0]=TempPoint.getx();
			this->xyz[k][1]=TempPoint.gety();
			this->xyz[k][2]=TempPoint.getz();
		}
	}

	BuildSigmaStruct(PRE_SIGMA_N,this->xyz);

	sort(PRE_SIGMA_N.begin(),PRE_SIGMA_N.end(),SIGMANCOMP);

	k=(int)PRE_SIGMA_N.size();
	m=-1;
	j=0;
	for(i=0;i<k;i++){
		if(m!=PRE_SIGMA_N[i].gtn){
			j++;
			m=PRE_SIGMA_N[i].gtn;
		}
	}
	SIGMA_N.resize(j);
	m=-1;
	j=0;
	for(i=0;i<k;i++){
		if(m!=PRE_SIGMA_N[i].gtn){
			SIGMA_N[j]=PRE_SIGMA_N[i];
			j++;
			m=PRE_SIGMA_N[i].gtn;
		}
	}

	this->n_nodes_dc = (int)SIGMA_N.size();
	this->n_nodes_c = this->n_nodes - this->n_nodes_dc;

	k = 0;	
	m = 0;
	j = this->n_nodes_c;
	for (i=0; i<kuzlov; i++)
	{
		if (isNodeInSubdomain[i])
		{
			if(m<this->n_nodes_dc && i==SIGMA_N[m].gtn){
				this->renumNodeFromNewToOld[j] = i;
				renumNodeFromOldToNew[i] = j;
				m++;
				j++;
			}
			else{
				this->renumNodeFromNewToOld[k] = i;
				renumNodeFromOldToNew[i] = k;
				k++;
			}
		}
	}

	m = this->n_nodes_dc;
	for (i=0; i<m; i++){
		SIGMA_N[i].gtn=renumNodeFromOldToNew[SIGMA_N[i].gtn];
		for(j=0;j<SIGMA_N[i].N;j++){
			SIGMA_N[i].bf[j]=renumNodeFromOldToNew[SIGMA_N[i].bf[j]];
		}
	}

	if (this->xyz_r)  {delete [] this->xyz_r;  this->xyz_r = NULL;}

	if (this->xyz)  {delete [] this->xyz;  this->xyz = NULL;}
	this->xyz = new double[this->n_nodes][3];

	for (i=0; i<this->n_nodes; i++)
	{
		pv::Point3D TempPoint;
		TempPoint=TaskCalcMesh->GetNode(this->renumNodeFromNewToOld[i]);
		this->xyz[i][0]=TempPoint.getx();
		this->xyz[i][1]=TempPoint.gety();
		this->xyz[i][2]=TempPoint.getz();
	}

	if (this->nver) {delete [] this->nver; this->nver = NULL;}
	this->nver = new long[this->n_elem][14];

	for (i=0; i<this->n_elem; i++){

		int elem = this->renumElemFromNewToOld[i];

		this->nver[i][13] = TaskCalcMesh->GetTypeOfElement(elem);

		for (j=0; j<8; j++){
			this->nver[i][j] = renumNodeFromOldToNew[TaskCalcMesh->GetNodeNumberOnElement(elem,j)];
		}
	}

	logfile<<"Resultant mesh: "<<this->n_nodes<<" nodes; "<<this->n_elem<<" elements; "<<this->n_nodes-this->n_nodes_c<<" terminal nodes;"<<'\n';
	cout<<"Resultant mesh: "<<this->n_nodes<<" nodes; "<<this->n_elem<<" elements; "<<this->n_nodes-this->n_nodes_c<<" terminal nodes;"<<'\n';

	int ktuzl,nc,p,ii,jj,kk,p1,p2;
	bool Correct_T_Flag;
	int pnt[2];
	ktuzl=this->n_nodes_dc;
	nc=this->n_nodes_c;
	do{
		Correct_T_Flag=false;
		for(i=0;i<ktuzl;i++)
		{
			double_eps a1,a2,a3;
			double_eps b1,b2,b3;

			if(SIGMA_N[i].N==2)
			{
				pnt[0]=SIGMA_N[i].bf[0];
				pnt[1]=SIGMA_N[i].bf[1];

				a1=this->xyz[pnt[0]][0];
				a2=this->xyz[pnt[0]][1];
				a3=this->xyz[pnt[0]][2];

				b1=this->xyz[pnt[1]][0];
				b2=this->xyz[pnt[1]][1];
				b3=this->xyz[pnt[1]][2];

				ii=(a1==b1);
				jj=(a2==b2);
				kk=(a3==b3);

				p = (ii && jj)? 2 : (ii && kk)? 1 : (jj && kk)? 0 : -1;
				p1 = (ii && jj) ? 0 : (ii && kk) ? 0 : (jj && kk) ? 1 : -1;
				p2 = (ii && jj) ? 1 : (ii && kk) ? 2 : (jj && kk) ? 2 : -1;
				if(p!=-1){
						Correct_T_Node(nc,SIGMA_N,i,p,pnt,Correct_T_Flag,this->xyz[pnt[0]][p1],this->xyz[pnt[0]][p2]);
					SIGMA_N[i].bf[0]=pnt[0];
					SIGMA_N[i].bf[1]=pnt[1];
					SIGMA_N[i].val[0]=(this->xyz[pnt[1]][p]-this->xyz[SIGMA_N[i].gtn][p])/(this->xyz[pnt[1]][p]-this->xyz[pnt[0]][p]);
					SIGMA_N[i].val[1]=(this->xyz[SIGMA_N[i].gtn][p]-this->xyz[pnt[0]][p])/(this->xyz[pnt[1]][p]-this->xyz[pnt[0]][p]);
				}
			}
		}
	}while(Correct_T_Flag);

	BuildTMatrix(SIGMA_N);

	int *p_jg_t;

	p_jg_t = (this->jg_t.size())? &(this->jg_t[0]) : NULL;

	this->p = new Portret(this->nver, this->n_elem, this->n_nodes, this->n_nodes_c,(long *) &(this->ig_t[0]),(long *) p_jg_t);
	this->p->Gen_T_Portrait2();

	vLV.resize(n_elem);

	for (i=0; i<this->n_nodes; i++)
	{
		pv::Point3D TempPoint;
		TempPoint=TaskCalcMesh->GetNodeTrue(this->renumNodeFromNewToOld[i]);
		this->xyz[i][0]=TempPoint.getx();
		this->xyz[i][1]=TempPoint.gety();
		this->xyz[i][2]=TempPoint.getz();
	}

	AsmGlobalMatrix();

	return 0;
}

void Subdomain::AsmGlobalMatrix()
{
	long i, j, k, m, it, jt, i_mu, j_nu;
	long ii, jj;
	long ig_n_1 = p->ig[n_nodes_c];

	if ((pr = new double[n_nodes_c*npls*ncmp]) == 0) Memory_allocation_error("pr", "AsmGlobalMatrix");
	if ((x = new double[n_nodes_c*npls*ncmp]) == 0) Memory_allocation_error("x", "AsmGlobalMatrix");
	if ((di = new double[n_nodes_c]) == 0) Memory_allocation_error("di", "AsmGlobalMatrix");
	if ((gg = new double[ig_n_1]) == 0) Memory_allocation_error("gg", "AsmGlobalMatrix");

	for (i=0; i<n_nodes_c; i++)
	{
		di[i] = 0;
	}

	for (i=0; i<n_nodes_c*npls*ntimes*ncmp; i++)
	{
		pr[i] = 0;
		x[i] = 0;
	}

	for (i=0; i<ig_n_1; i++)
		gg[i] = 0;

	for(i=0; i<n_elem; i++)
	{
		Hex_Local_Matrix L(i, nver, xyz);
		L.CalcMassMatrix();

		for(j=0; j<8; j++)
		{
			vLV[i].g[j]=0.0;
			for(k=0; k<8; k++){vLV[i].g[j]+=L.b[j][k];}

			ii = nver[i][j];

			if(ii < n_nodes_c)
			{ 
				di[ii] += L.b[j][j];
			}
			else
			{
				for(it = p->ig_t[ii]; it <= p->ig_t[ii+1]-1; it++)
				{
					i_mu = p->jg_t[it];

					for(jt = p->ig_t[ii]; jt <= p->ig_t[ii+1]-1; jt++)
					{
						j_nu = p->jg_t[jt];

						if(j_nu < i_mu)
						{
							for(m = p->ig[i_mu]; m <= p->ig[i_mu+1]-1; m++)
							{
								if(p->jg[m] == j_nu)
									gg[m] += L.b[j][j] * gg_t[it] * gg_t[jt];
							}
						}
						else if(j_nu == i_mu)
						{
							di[i_mu] += L.b[j][j] * gg_t[it] * gg_t[jt];
						}
					}
				}
			}

			for(k=0; k<8; k++)
			{
				jj = nver[i][k];

				if(ii < n_nodes_c && jj < n_nodes_c)
				{ 
					if(jj < ii) 
					{
						for(m=p->ig[ii]; m<=p->ig[ii+1]-1; m++)
						{
							if(p->jg[m] == jj)
								gg[m] += L.b[j][k];
						}
					}
				}

				else if(ii >= n_nodes_c && jj < n_nodes_c && j != k)
				{
					for(it = p->ig_t[ii]; it <= p->ig_t[ii+1]-1; it++)
					{
						i_mu = p->jg_t[it];

						if(jj < i_mu)
						{
							for(m=p->ig[i_mu]; m<=p->ig[i_mu+1]-1; m++)
								if(p->jg[m] == jj)
									gg[m] += L.b[j][k] * gg_t[it];
						}
						else if(jj == i_mu)
						{
							di[i_mu] += L.b[j][k] * gg_t[it];
						}
					}
				}
				else if(ii < n_nodes_c && jj >= n_nodes_c && j!=k)
				{
					for(it = p->ig_t[jj]; it <= p->ig_t[jj+1]-1; it++)
					{
						j_nu = p->jg_t[it];
						if(j_nu < ii)
						{
							for(m = p->ig[ii]; m <= p->ig[ii+1]-1; m++)
								if(p->jg[m] == j_nu)
									gg[m] += L.b[j][k] * gg_t[it];
						}
						else if(j_nu == ii)
						{
							di[j_nu] += L.b[j][k] * gg_t[it];
						}
					}
				}
				else if(ii >= n_nodes_c && jj >= n_nodes_c /*&& ii!=jj*/ && j!=k)
				{
					for(it = p->ig_t[ii]; it <= p->ig_t[ii+1]-1; it++)
					{
						i_mu = p->jg_t[it];
						for(jt = p->ig_t[jj]; jt <= p->ig_t[jj+1]-1; jt++)
						{
							j_nu = p->jg_t[jt];
							if(j_nu < i_mu)
							{
								for(m = p->ig[i_mu]; m <= p->ig[i_mu+1]-1; m++)
									if(p->jg[m] == j_nu)
										gg[m] += L.b[j][k] * gg_t[it] * gg_t[jt];
							}
							else if(j_nu == i_mu)
							{
								di[i_mu] += L.b[j][k] * gg_t[it] * gg_t[jt];
							}
						}
					}
				}
			}
		}
	}

	prds.factorize_rsf(n_nodes_c,(int *)p->ig,(int *)p->jg,gg,di,nthreads);
}

void Subdomain::AsmGlobalVector(int ipls,int icmp)
{
	long i, j, it, i_mu;
	long ii;
	long ig_n_1 = p->ig[n_nodes_c];

	for (i=0; i<n_nodes_c; i++)
		pr[(icmp*npls+ipls)*n_nodes_c+i] = 0;

	for(i=0; i<n_elem; i++)
	{
		Local_Vector &lv=vLV[i];

		for(j=0; j<8; j++)
		{
			ii = nver[i][j];

			if(ii < n_nodes_c)
			{
				pr[(icmp*npls+ipls)*n_nodes_c+ii] += lv.g[j]*ValueInCenter[i];
			}
			else
			{
				for(it = p->ig_t[ii]; it<=p->ig_t[ii+1]-1; it++)
				{
					i_mu = p->jg_t[it];
					pr[(icmp*npls+ipls)*n_nodes_c+i_mu] += lv.g[j]*ValueInCenter[i]*gg_t[it];
				}
			}
		}
	}
}
void Subdomain::CalcRightPartVect(T_Brick &L)
{
	int i,j;
	double BzOnElem;

	BzOnElem=ValueInCenter[L.num];

	Hex_Local_Matrix LM(L.num,nver,xyz);
	if(LM.type_of_hex<31){
		LM.Calc_local_matrix_b_for_parallelepiped();
	}
	else{
		int i1, j1;
		double gauss_3_mult;

		for(i=0; i<8; i++)
		{
			for(j=0; j<8; j++)
				LM.b[i][j] = 0.0;
		}

		for(i=0; i<27; i++)
		{
			LM.Calc_J(i);
			gauss_3_mult = gauss_3_A_all[i]*LM.det_J_abs;

			for(i1=0; i1<8; i1++)
			{
				for(j1=0; j1<8; j1++)
					LM.b[i1][j1] += gauss_3_phi[i][i1]*gauss_3_phi[i][j1]*gauss_3_mult;
			}
		}
	}

	for(i=0;i<8;i++){
		L.g8[i]=0;
		for(j=0;j<8;j++){
			L.g8[i]+=LM.b[i][j];
		}
		L.g8[i]*=BzOnElem;
	}
}

void Subdomain::build_elem_neib_first(int nx,int ny,int nz,int *regular,vector<bool> &ElemCheck)
{
	int i,j,k,n_el,m_el,p_el,t_el;
	int ex,ey,ez;
	int bx,by,bz;
	int m,size;

	ElemNeibVec.clear();
	ElemNeibVec.resize(n_elem);

	ex=nx-1;ey=ny-1;ez=nz-1;
	bx=ex-1;by=ey-1;bz=ez-1;

	for(k=0;k<ez;k++){
		for(j=0;j<ey;j++){
			for(i=0;i<ex;i++){

				n_el=regular[((ey*k)+j)*ex+i];
				p_el=renumElemFromOldToNew[n_el];

				if(ElemCheck[n_el]){

					if(i>0){
						m_el=regular[((ey*k)+j)*ex+(i-1)];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[0].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[0][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[0].push_back(t_el);
						}
					}

					if(j>0){
						m_el=regular[((ey*k)+(j-1))*ex+i];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[1].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[1][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[1].push_back(t_el);
						}
					}

					if(k>0){
						m_el=regular[((ey*(k-1))+j)*ex+i];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[2].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[2][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[2].push_back(t_el);
						}						
					}

					if(i<bx){
						m_el=regular[((ey*k)+j)*ex+(i+1)];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[3].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[3][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[3].push_back(t_el);
						}
					}

					if(j<by){
						m_el=regular[((ey*k)+(j+1))*ex+i];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[4].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[4][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[4].push_back(t_el);
						}
					}

					if(k<bz){
						m_el=regular[((ey*(k+1))+j)*ex+i];
						t_el=renumElemFromOldToNew[m_el];
						if(n_el!=m_el){
							size=(int)ElemNeibVec[p_el].neib[5].size();
							for(m=0;m<size;m++){
								if(t_el==ElemNeibVec[p_el].neib[5][m])break;
							}
							if(m==size && ElemCheck[m_el])ElemNeibVec[p_el].neib[5].push_back(t_el);
						}						
					}
				}
			}
		}
	}
}

void Subdomain::BuildSigmaStruct(vector<_SIGMA_N> &SIGMA_N,double (*xyz_p)[3])
{
	int i,ii,jj,ll,nn,mm,rr,tt,p,p1,p2;
	_SIGMA_N TTNE;

	for(i=0;i<n_elem;i++){
		int node_new[4],node_old[4];
		double_eps a1,b1;
		double_eps c1,c2,d1,d2;
		double invhc,invhd,invhchd;
		int eq;
		nn=renumElemFromNewToOld[i];
		for(ii=0;ii<6;ii++){
			switch(ii){
				case 0 : {ll=3;p=0;p1=1;p2=2;break;}
				case 1 : {ll=4;p=1;p1=0;p2=2;break;}
				case 2 : {ll=5;p=2;p1=0;p2=1;break;}
				case 3 : {ll=0;p=0;p1=1;p2=2;break;}
				case 4 : {ll=1;p=1;p1=0;p2=2;break;}
				case 5 : {ll=2;p=2;p1=0;p2=1;break;}
			}
			node_new[0]=nver[nn][Sn[ii][0]];
			node_new[1]=nver[nn][Sn[ii][1]];
			node_new[2]=nver[nn][Sn[ii][2]];
			node_new[3]=nver[nn][Sn[ii][3]];

			tt=(int)ElemNeibVec[i].neib[ii].size();
			for(rr=0;rr<tt;rr++){
				mm=renumElemFromNewToOld[ElemNeibVec[i].neib[ii][rr]];

				node_old[0]=nver[mm][Sn[ll][0]];
				node_old[1]=nver[mm][Sn[ll][1]];
				node_old[2]=nver[mm][Sn[ll][2]];
				node_old[3]=nver[mm][Sn[ll][3]];

				c1=xyz_p[node_old[0]][p1];
				c2=xyz_p[node_old[3]][p1];
				d1=xyz_p[node_old[0]][p2];
				d2=xyz_p[node_old[3]][p2];
				invhc=1.0/((c2-c1).get_value());
				invhd=1.0/((d2-d1).get_value());
				invhchd=invhc*invhd;

				for(jj=0;jj<4;jj++){
					a1=xyz_p[node_new[jj]][p1];
					b1=xyz_p[node_new[jj]][p2];
					TTNE.gtn=node_new[jj];
					eq=(a1==c1 || a1==c2)+(b1==d1 || b1==d2);
					if(eq!=2 && a1>=c1 && a1<=c2 && b1>=d1 && b1<=d2){
						if(b1==d1){
							TTNE.N=2;
							TTNE.bf[0]=node_old[0];
							TTNE.bf[1]=node_old[1];
							TTNE.val[0]=(c2-a1).value*invhc;
							TTNE.val[1]=(a1-c1).value*invhc;
							SIGMA_N.push_back(TTNE);
						}
						else if(b1==d2){
							TTNE.N=2;
							TTNE.bf[0]=node_old[2];
							TTNE.bf[1]=node_old[3];
							TTNE.val[0]=(c2-a1).value*invhc;
							TTNE.val[1]=(a1-c1).value*invhc;
							SIGMA_N.push_back(TTNE);
						}
						else if(a1==c1){
							TTNE.N=2;
							TTNE.bf[0]=node_old[0];
							TTNE.bf[1]=node_old[2];
							TTNE.val[0]=(d2-b1).value*invhd;
							TTNE.val[1]=(b1-d1).value*invhd;
							SIGMA_N.push_back(TTNE);
						}
						else if(a1==c2){
							TTNE.N=2;
							TTNE.bf[0]=node_old[1];
							TTNE.bf[1]=node_old[3];
							TTNE.val[0]=(d2-b1).value*invhd;
							TTNE.val[1]=(b1-d1).value*invhd;
							SIGMA_N.push_back(TTNE);
						}
						else{
							TTNE.N=4;
							TTNE.bf[0]=node_old[0];
							TTNE.bf[1]=node_old[1];
							TTNE.bf[2]=node_old[2];
							TTNE.bf[3]=node_old[3];
							TTNE.val[0]=(c2-a1).value*(d2-b1).value*invhchd;
							TTNE.val[1]=(a1-c1).value*(d2-b1).value*invhchd;
							TTNE.val[2]=(c2-a1).value*(b1-d1).value*invhchd;
							TTNE.val[3]=(a1-c1).value*(b1-d1).value*invhchd;
							SIGMA_N.push_back(TTNE);
						}
					}
				}	
			}
		}
	}
}

void Subdomain::Correct_T_Node(int nc,vector<_SIGMA_N> &SIGMA_N,int i,int p,int *pnt,bool &Correct_T_Flag, double pVal1, double pVal2)
{
	int j,k,m,n;
	int ii,jj,kk,t,t1,t2;
	int tnt[2];
	n=nc+i;
	
	for(j=0;j<SIGMA_N[i].N;j++)
	{
		double_eps a1,a2,a3;
		double_eps b1,b2,b3;

		k=SIGMA_N[i].bf[j];
		m=k-nc;

		if(k>=nc && SIGMA_N[m].N==2)
		{
			tnt[0]=SIGMA_N[m].bf[0];
			tnt[1]=SIGMA_N[m].bf[1];
	
			a1=xyz[tnt[0]][0];
			a2=xyz[tnt[0]][1];
			a3=xyz[tnt[0]][2];

			b1=xyz[tnt[1]][0];
			b2=xyz[tnt[1]][1];
			b3=xyz[tnt[1]][2];

			ii=(a1==b1);
			jj=(a2==b2);
			kk=(a3==b3);

			t = (ii && jj)? 2 : (ii && kk)? 1 : (jj && kk)? 0 : -1;
			t1 = (ii && jj) ? 0 : (ii && kk) ? 0 : (jj && kk) ? 1 : -1;
			t2 = (ii && jj) ? 1 : (ii && kk) ? 2 : (jj && kk) ? 2 : -1;

			if(p==t && fabs(xyz[tnt[0]][t1]-pVal1)<EpsComp && fabs(xyz[tnt[0]][t2]-pVal2)<EpsComp)
			{
				if(xyz[tnt[1]][p]<xyz[tnt[0]][p])
				{
					int tmpi = tnt[1];
					tnt[1] = tnt[0];
					tnt[0] = tmpi;
				}

				if(xyz[tnt[0]][p]<xyz[pnt[0]][p]){
					pnt[0]=tnt[0];
					if(pnt[0]>=nc && SIGMA_N[pnt[0]-nc].N==2)Correct_T_Node(nc,SIGMA_N,pnt[0]-nc,p,pnt,Correct_T_Flag,pVal1,pVal2);
					Correct_T_Flag=true;
				}
				if(xyz[tnt[1]][p]>xyz[pnt[1]][p]){
					pnt[1]=tnt[1];
					if(pnt[1]>=nc && SIGMA_N[pnt[1]-nc].N==2)Correct_T_Node(nc,SIGMA_N,pnt[1]-nc,p,pnt,Correct_T_Flag,pVal1,pVal2);
					Correct_T_Flag=true;
				}
			}
		}
	}
}

void Subdomain::BuildTMatrix(vector<_SIGMA_N> &SIGMA_N)
{
	int a,b,ci;
	double cd;
	int j,li,t,m;
	double vklad;
	vector<double> dv;
	vector<int> iv;

	ig_t.clear();
	jg_t.clear();

	m=0;
	ig_t.resize(n_nodes_c+1,m);
	for(j=0;j<n_nodes_dc;j++){
		t=0;
		for(li=0;li<SIGMA_N[j].N;li++){
			vklad=SIGMA_N[j].val[li];
			calc_Tij(SIGMA_N,j,li,iv,dv,t,vklad);
		}
		t=(int)iv.size();

		for(a=0;a<t;a++){
			gg_t.push_back(dv[a]);
			jg_t.push_back(iv[a]);
			m++;
		}
		ig_t.push_back(m);
		iv.clear();
		dv.clear();
	}
}

void Subdomain::calc_Tij(vector<_SIGMA_N> &SIGMA_N,int j,int li,vector<int> &iv,vector<double> &dv,int &t,double vklad)
{
	int i;
	i=SIGMA_N[j].bf[li];
	if(i<n_nodes_c){
		iv.push_back(i);
		dv.push_back(vklad);
		t++;
	}
	else{
		int lj;
		double vkl;
		for(lj=0;lj<SIGMA_N[i-n_nodes_c].N;lj++){
			vkl=vklad*SIGMA_N[i-n_nodes_c].val[lj];
			calc_Tij(SIGMA_N,i-n_nodes_c,lj,iv,dv,t,vkl);
		}
	}
}

void Subdomain::CalcValuesAll_Node(double *v)
{
	long i, j;
	double s;

	for (i=n_nodes_c; i<n_nodes; i++)
	{
		s = 0.0;
		for (j=ig_t[i]; j<=ig_t[i+1]-1; j++)
			s += v[jg_t[j]]*gg_t[j];
		v[i] = s;
	}
}
