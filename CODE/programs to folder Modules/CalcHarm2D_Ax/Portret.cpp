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
 *  This file contains code of some basic routines for construction of a matrix portrait of a finite element SLAE 
 * 
 * 
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin   
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                      
 *  Version 2.0 December 10, 2024                                                         
*/ 

#include "stdafx.h"
#include "Portret.h"
#include "in_out.h"
#include "For_Solvers.h"

void Portret::Add_In_Ordered_List(list *s, long x)
{
	list *head, *cur, *prev, *list_new;

	head = s->next;
	cur = head;
	prev = NULL;

	while(cur != NULL)
	{
   		if(x < cur->number)
      		break;
		prev = cur;
		cur = cur->next;
	}

	if(prev == NULL)
	{
	  	list_new = new list;
		if(list_new == 0)
		{
			char var[] = {"list_new"};
			char func[] = {"Portret::Add_In_Ordered_List"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		list_new->number = x;

		list_new->next = head;
		s->next = list_new;
	}
	else if(prev->number!=x)
	{
	  	list_new = new list;
		if(list_new == 0)
		{
			char var[] = {"list_new"};
			char func[] = {"Portret::Add_In_Ordered_List"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		list_new->number = x;

   		list_new->next = prev->next;
		prev->next = list_new;
	}
}

Portret::Portret(long (*nvtr)[4], long n_of_elements, long n_of_nodes, long n_c)
{
	this->nvtr4 = nvtr;
	this->n_of_elements = n_of_elements;
	this->is_mem_ig_jg_allocated = false;
	this->n = n_of_nodes;
	this->n_c = n_c;
}

Portret::Portret(long (*nvtr)[8], long n_of_elements, long n_of_nodes, long n_c)
{
	this->nvtr = nvtr;
	this->n_of_elements = n_of_elements;
	this->is_mem_ig_jg_allocated = false;
	this->n = n_of_nodes;
	this->n_c = n_c;
}

Portret::Portret(long (*nvtr4)[4], long n_of_elements, long n_of_nodes)
{
	this->nvtr4 = nvtr4;
	this->n_of_elements = n_of_elements;
	this->is_mem_ig_jg_allocated = false;
	this->n = n_of_nodes;
}

Portret::Portret(long (*nver)[14], long n_of_elements, long n_of_nodes, long n_c)
{
	this->nver = nver;
	this->n_of_elements = n_of_elements;
	this->is_mem_ig_jg_allocated = false;
	this->n = n_of_nodes;
	this->n_c = n_c;
}

Portret::Portret()
{
	this->is_mem_ig_jg_allocated = false;
}

void Portret::Clear_List(list *s)
{
	list *cur_pos, *next_pos;

	cur_pos = s;

	while(cur_pos != NULL)
	{
		next_pos = cur_pos->next;
		delete cur_pos;
		cur_pos = next_pos;
	}
}

Portret::~Portret()
{
	if(this->is_mem_ig_jg_allocated==true)
	{
		delete [] this->ig;
		delete [] this->jg;
	}
}

void Portret::Gen_T_Portrait()
{
	Portret P;
	long i, j, k;
	long tmp1, tmp2, e;
	list *l, *s;
	std::vector<long> local, loc;
	In_Out R;

	printf("T_Portrait... ");

	s = new list[n_c];
	if(s == 0)
	{
		char var[] = {"s"};
		char func[] = {"Portret::Gen_T_Portrait"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<this->n_of_elements; i++)
	{
		for(j=0; j<8; j++)
		{
			e = nvtr[i][j];
			local.push_back(e);
		}

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		for(j=0; j<(long)loc.size(); j++)
		for(k=0; k<(long)loc.size(); k++)
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
		tmp1 = (long)loc.size();
		for(j=0; j<tmp1; j++) loc.pop_back();
	}

	size_jg = 0;
	for(i=0; i<n_c; i++)
		size_jg += P.Size_Of_List(s[i].next);

	ig = new long[n_c + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portret::Gen_T_Portrait"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	jg = new long[size_jg];
	if(jg == 0)
	{
		char var[] = {"jg"};
		char func[] = {"Portret::Gen_T_Portrait"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	is_mem_ig_jg_allocated = true;

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
	delete [] s;

	printf("done.\n");
}

void Portret::Gen_T_Portrait2()
{
	Portret P;
	long i, j, k;
	long tmp1, tmp2, e;
	list *l, *s;
	std::vector<long> local, loc;
	In_Out R;

	if ((s = new list[n_c]) == 0) Memory_allocation_error("s", "Gen_T_Portrait2");

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<n_of_elements; i++)
	{
		for(j=0; j<8; j++)
		{
			e = nver[i][j];
			local.push_back(e);
		}

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		for(j=0; j<(long)loc.size(); j++)
		{
			for(k=0; k<(long)loc.size(); k++)
			{
				tmp1 = loc[k];
				tmp2 = loc[j];

				if(tmp1 < tmp2)
					P.Add_In_Ordered_List(&s[tmp2], tmp1);					
			}
		}

		loc.clear();
		local.clear();
	}

	size_jg = 0;
	for(i=0; i<n_c; i++)
		size_jg += P.Size_Of_List(s[i].next);

	if ((ig = new long[n_c + 1]) == 0) Memory_allocation_error("ig", "Gen_T_Portrait2");
	if ((jg = new long[size_jg]) == 0) Memory_allocation_error("jg", "Gen_T_Portrait2");;

	is_mem_ig_jg_allocated = true;

	i = 0;
	ig[0] = 0;
	for(k=0; k<n_c; k++)
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

	delete [] s;
}

void Portret::Gen_Portret_2d_Linear()
{
	long i, j, k;
	long tmp1, tmp2;
	list *l;

	printf("Portret_2d_Linear... ");

	s = new list[n];
	if(s == 0)
	{
		char var[] = {"s"};
		char func[] = {"Portret::Gen_Portret_2d_Linear"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	for(i=0; i<n; i++)
		s[i].next = NULL;

	for(i=0; i<this->n_of_elements; i++)
		for(j=0; j<4; j++)
			for(k=0; k<4; k++)
			{
				tmp1 = nvtr4[i][k];
				tmp2 = nvtr4[i][j];
				if(tmp1 < tmp2)
					this->Add_In_Ordered_List(&s[tmp2], tmp1);					
			}

		this->size_jg = 0;
		for(i=0; i<this->n; i++)
			this->size_jg += this->Size_Of_List(s[i].next);

		ig = new long[n+1];
		if(ig == 0)
		{
			char var[] = {"ig"};
			char func[] = {"Portret::Gen_Portret_2d_Linear"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		jg = new long[size_jg];
		if(jg == 0)
		{
			char var[] = {"jg"};
			char func[] = {"Portret::Gen_Portret_2d_Linear"};
			printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
			exit(1);
		}

		is_mem_ig_jg_allocated = true;

		i = 0;
		ig[0] = 0;
		for(k=0;k<this->n;k++)
		{
			l = s[k].next;
			while(l != NULL)
			{
				jg[i] = l->number;
				i++;
				l = l->next;
			}
			this->Clear_List(this->s[k].next);
			ig[k+1] = i;
		}
		delete [] s;

		printf("done.\n");
}

long Portret::Size_Of_List(list *s)
{
	long i;
	list *cur_pos;

	i = 0;
	cur_pos = s;

	while(cur_pos != NULL)
	{
		i++;
		cur_pos = cur_pos->next;
	}

	return i;
}
//------------------------------------------------------ 
// Constructing a porter of matrix of a finite element SLAE 
//------------------------------------------------------ 
void Portret::Gen_Portret_2d_rect_T()
{
	long i, j, k;
	long tmp1, tmp2, e;
	list *l, *s;
	std::vector<long> local, loc;
	In_Out R;

	printf("T_Portrait... ");

	s = new list[n_c];
	if(s == 0)
	{
		char var[] = {"s"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	for(i=0; i<n_c; i++)
		s[i].next = NULL;

	for(i=0; i<this->n_of_elements; i++)
	{
		for(j=0; j<4; j++)
		{
			e = nvtr4[i][j];
			local.push_back(e);
		}

		std::sort(local.begin(),local.end());
		std::unique_copy(local.begin(), local.end(), back_inserter(loc));

		for(j=0; j<(long)loc.size(); j++)
		{
			for(k=0; k<(long)loc.size(); k++)
			{
				tmp1 = loc[k];
				tmp2 = loc[j];
				if(tmp1 < tmp2)
				{
					Add_In_Ordered_List(&s[tmp2], tmp1);					
				}
			}
		}

		tmp1 = (long)local.size();
		for (j = 0; j < tmp1; j++) local.pop_back();
		tmp1 = (long)loc.size();
		for (j = 0; j < tmp1; j++) loc.pop_back();
	}

	size_jg = 0;
	for(i=0; i<n_c; i++)
		size_jg += Size_Of_List(s[i].next);

	ig = new long[n_c + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	jg = new long[size_jg];
	if(jg == 0)
	{
		char var[] = {"jg"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	is_mem_ig_jg_allocated = true;

	i = 0;
	ig[0] = 0;
	for(k=0;k<n_c;k++)
	{
		l = s[k].next;
		while(l != NULL)
		{
			jg[i] = l->number;
			i++;
			l = l->next;
		}
		Clear_List(s[k].next);
		ig[k+1] = i;
	}
	delete [] s;

	printf("done.\n");
}
