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
 *  Written by Prof. Yuri G. Soloveichik, Prof. Marina G. Persova, D.Sc. Denis V. Vagin   
 *  Novosibirsk State Technical University,                                               
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                     
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                      
 *  Version 2.0 December 10, 2024                                                         
 *
*/

#pragma once

struct list
{
	int number;
	list *next;
};

class Portret
{
public:

	int *ig;
	int *jg;
	list *s;
	int n;
	int n_of_elements;
	int size_jg;
	int (*nvtr4)[4];
	int n_c;
	
	Portret(int (*nvtr)[4], int n_of_elements, int n_of_nodes, int n_c);
	~Portret();

	void Gen_Portret_2d_rect_T_3x3();
	// Construction of the matrix port of the finite element SLAE
	void Gen_Portret_2d_rect_T();

	void Add_In_Ordered_List(list *s, int x);
	void Clear_List(list *s);
	int Size_Of_List(list *s);
};
