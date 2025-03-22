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
 *  This file contains headers for some basic routines for construction of a matrix portrait of a finite element SLAE
 *
 *
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin   
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                      
 *  Version 2.0 December 10, 2024                                                         
*/

#pragma once

struct list
{
	long number;
	list *next;
};

class Portret
{
public:

	long *ig;
	long *jg;
	list *s;
	bool is_mem_ig_jg_allocated;
	long n;
	long n_of_elements;
	long size_jg;
	long (*nvtr)[8];
	long (*nvtr4)[4];
	long (*nver)[14];
	long n_c;
	Portret(long (*nver)[14], long n_of_elements, long n_of_nodes, long n_c);
	Portret(long (*nvtr)[4],  long n_of_elements, long n_of_nodes, long n_c);
	Portret(long (*nvtr)[8],  long n_of_elements, long n_of_nodes, long n_c);
	Portret(long (*nvtr4)[4], long n_of_elements, long n_of_nodes);
	Portret();
	~Portret();
	void Gen_T_Portrait();
	void Gen_T_Portrait2();
	void Gen_Portret_2d_Linear();
	// Constructing a porter of matrix of a finite element SLAE 
	void Gen_Portret_2d_rect_T();
	void Add_In_Ordered_List(list *s, long x);
	void Clear_List(list *s);
	long Size_Of_List(list *s);
};
