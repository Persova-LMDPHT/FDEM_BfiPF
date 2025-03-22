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
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin   
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                      
 *  Version 2.0 December 10, 2024                                                         
 *  
*/

#include "stdafx.h"
#include "Portret.h"
#include "in_out.h"

void Portret::Add_In_Ordered_List(list *s, int x)
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

Portret::Portret(int (*nvtr)[4], int n_of_elements, int n_of_nodes, int n_c)
{
	this->nvtr4 = nvtr;
	this->n_of_elements = n_of_elements;
	this->n = n_of_nodes;
	this->n_c = n_c;
	ig = NULL;
	jg = NULL;
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
	if (ig) {delete [] ig; ig = NULL;}
	if (jg) {delete [] jg; jg = NULL;}
}

int Portret::Size_Of_List(list *s)
{
	int i;
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

void Portret::Gen_Portret_2d_rect_T_3x3()
{
	int i, j, k, j1, k1;
	int tmp1, tmp2, e;
	list *l, *s;
	std::vector<int> local, loc;
	In_Out R;

	printf("T_Portrait... ");

	s = new list[n_c*3];
	if(s == 0)
	{
		char var[] = {"s"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	for(i=0; i<n_c*3; i++)
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

		for(j=0; j<(int)loc.size(); j++)
			for(k=0; k<(int)loc.size(); k++)
				for(j1=0; j1<3; j1++)
					for(k1=0; k1<3; k1++)
					{
						tmp1 = loc[k]*3 + k1;
						tmp2 = loc[j]*3 + j1;
						if(tmp1 < tmp2)
						{
							Add_In_Ordered_List(&s[tmp2], tmp1);					
						}
					}

					tmp1 = (int)local.size();
					for(j=0; j<tmp1; j++) local.pop_back();
					tmp1 = (int)loc.size();
					for(j=0; j<tmp1; j++) loc.pop_back();
	}

	size_jg = 0;
	for(i=0; i<n_c*3; i++)
		size_jg += Size_Of_List(s[i].next);

	ig = new int[n_c*3 + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	jg = new int[size_jg];
	if(jg == 0)
	{
		char var[] = {"jg"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	i = 0;
	ig[0] = 0;
	for(k=0;k<n_c*3;k++)
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
//------------------------------------------------------ 
// Construction of the matrix port of the finite element SLAE
//------------------------------------------------------ 
void Portret::Gen_Portret_2d_rect_T()
{
	int i, j, k;
	int tmp1, tmp2, e;
	list *l, *s;
	std::vector<int> local, loc;
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

		for(j=0; j<(int)loc.size(); j++)
		{
			for(k=0; k<(int)loc.size(); k++)
			{
				tmp1 = loc[k];
				tmp2 = loc[j];
				if(tmp1 < tmp2)
				{
					Add_In_Ordered_List(&s[tmp2], tmp1);					
				}
			}
		}

		tmp1 = (int)local.size();
		for(j=0; j<tmp1; j++) local.pop_back();
		tmp1 = (int)loc.size();
		for(j=0; j<tmp1; j++) loc.pop_back();
	}

	size_jg = 0;
	for(i=0; i<n_c; i++)
		size_jg += Size_Of_List(s[i].next);

	ig = new int[n_c + 1];
	if(ig == 0)
	{
		char var[] = {"ig"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

	jg = new int[size_jg];
	if(jg == 0)
	{
		char var[] = {"jg"};
		char func[] = {"Portret::Gen_Portret_2d_rect_T"};
		printf("Memory allocation error for variable \"%s\" in function \"%s\"\n",var,func);
		exit(1);
	}

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
