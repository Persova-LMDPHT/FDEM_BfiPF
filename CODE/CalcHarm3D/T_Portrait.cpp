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
 *  This file contains code of functions for construction of a matrix portrait of 3D VFEM SLAE
 * 
 * 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/T_Portrait.cpp
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)
 *  
*/

#include "stdafx.h"  
#include "T_Portrait.h"
#include "Portret.h"
extern ofstream logfile;

T_Portrait::T_Portrait(long *ed, long n, long n_c, long n_elem)
{
	this->n = n;
	this->n_c = n_c;
	this->n_elem = n_elem;
	this->ed = (long(*)[25])ed;

	ig = NULL;
	jg = NULL;
	idi = NULL;
	ijg = NULL;
}

T_Portrait::~T_Portrait()
{
	if(ig) {delete [] ig; ig=NULL;}
	if(jg) {delete [] jg; jg=NULL;}
	if(idi) { delete [] idi; idi=NULL; }
	if(ijg) { delete [] ijg; ijg=NULL; }
}

void T_Portrait::Gen_Portrait()
{
	Portret P;
	long i, j, k;
	long tmp1, tmp2, e;
	long loc_size;
	list *l=NULL, *s=NULL;
	std::vector<long> local, loc;

	logfile << "T_Portrait... ";

	s = new list[n_c];
	if(s==NULL)
		Memory_allocation_error("s", "T_Portrait::Gen_Portrait");

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<this->n_elem; i++)
	{
		for(j=0; j<12; j++)
		{
			e = ed[i][j];
			local.push_back(e);
		}

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		loc_size = (long)loc.size();

		for(j=0; j<loc_size; j++)
		for(k=0; k<loc_size; k++)
		{
			tmp1 = loc[k];
			tmp2 = loc[j];
			if(tmp1 < tmp2)
			{
				P.Add_In_Ordered_List(&s[tmp2], tmp1);
			}
		}

		tmp1 = (long)local.size();
		for(j=0; j<tmp1; j++) local.pop_back();
		for(j=0; j<loc_size; j++) loc.pop_back();
	}

	this->size_jg = 0;
	for(i=0; i<this->n_c; i++)
		this->size_jg += P.Size_Of_List(s[i].next);

	ig = new long[this->n_c + 1];
	if(ig == 0)
		Memory_allocation_error("ig", "T_Portrait::Gen_Portrait");

	jg = new long[this->size_jg];
	if(jg == 0)
		Memory_allocation_error("jg", "T_Portrait::Gen_Portrait");

	i = 0;
	ig[0] = 0;
	for(k=0;k<this->n_c;k++)
	{
		l = s[k].next;
		while(l != NULL)
		{
			jg[i] = l->number;
			i++;
			l = l->next;
		}
		P.Clear_List(s[k].next);
		ig[k+1] = i;
	}

	if(s) {delete [] s; s=NULL;}

	logfile << "done.\n";
}

void T_Portrait::Gen_idi_ijg(long *nvkat, long (*nver)[14])
{
	long i, j, k, m, it, jt, ii, jj, i_mu, j_nu;
	long i2, j2;
	long t, t2;
	long a[12][12];

	if((idi = new long[n_c+1])==0) Memory_allocation_error("idi", "T_Portrait::Gen_idi_ijg");
	if((ijg = new long[size_jg+1])==0) Memory_allocation_error("ijg", "T_Portrait::Gen_idi_ijg");

	for(i=0; i<n_c+1; i++)
		idi[i] = 0;

	for(i=0; i<size_jg+1; i++)
		ijg[i] = 0;

	for (i=0; i<n_elem; i++)
	{
		if(false)
		{
			for(i2=0; i2<12; i2++)
				for(j2=0; j2<12; j2++)
					a[i2][j2] = 1;
		}
		else
		{
			if (nver[i][13] > 30)
			{
				for(i2=0; i2<12; i2++)
					for(j2=0; j2<12; j2++)
						a[i2][j2] = 2;
			} 
			else
			{
				for(i2=0; i2<12; i2++)
					for(j2=0; j2<12; j2++)
						a[i2][j2] = FILTER_MASS_MATRIX_VEC[i2][j2];
			}
		}

		for(j=0; j<12; j++)
		{
			ii = ed[i][j];

			if(ii < n_c)
				Set_type_of_block(idi, ii, a[j][j]);

			for(k=0; k<12; k++)
			{
				jj = ed[i][k];

				if(jj < ii)
				{
					for(m=ig[ii]; m<=ig[ii+1]-1; m++)
					{
						if(jg[m]==jj)
						{
							Set_type_of_block(ijg, m, a[j][k]);
							break;
						}
					}
				}
			}
		}
	}

	t = idi[0];
	idi[0] = 0;

	for (i=0; i<n_c; i++)
	{
		t2 = idi[i] + t;
		t = idi[i+1];
		idi[i+1] = t2;
	}

	t = ijg[0];
	ijg[0] = 0;

	for (i=0; i<size_jg; i++)
	{
		t2 = ijg[i] + t;
		t = ijg[i+1];
		ijg[i+1] = t2;
	}
}

void T_Portrait::Set_type_of_block(long *target_block, long adr, long type)
{
	if(target_block[adr]!=2)
		target_block[adr] = type;
}
