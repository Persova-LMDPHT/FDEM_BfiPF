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
 *  This file contains headers of functions for construction of a matrix portrait of 3D VFEM SLAE
 * 
 * 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/BuildMatrix/T_Portrait.h          
 *  Novosibirsk State Technical University,                                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                               
 *  
*/

#pragma once
class T_Portrait
{
public:
	long (*ed)[25];
	long n_elem;
	long n;
	long n_c;
	long *ig; 
	long *jg;
	long size_jg;
	long *idi;
	long *ijg;
	T_Portrait(long *ed, long n, long n_c, long n_elem);
	~T_Portrait();
	void Gen_Portrait();
	void Gen_idi_ijg(long *nvkat, long (*nver)[14]);
	void Set_type_of_block(long *target_array, long adr, long type);
};
const int FILTER_MASS_MATRIX_VEC[12][12] = { 
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,  
   1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2
};
