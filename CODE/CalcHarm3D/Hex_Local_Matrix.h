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
 *  This file contains headers for computation of local matrices on hexahedrons
 * 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/Hex_Local_Matrix.h      
 *  Novosibirsk State Technical University,                                                                      
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                            
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                             
 *  
*/  

#pragma once

class Hex_Local_Matrix
{
public:
	int number_of_element;
	double x[8], y[8], z[8];   // coordinates of the vertices of the hexahedrons	    
	double J[3][3];            // the Jacobi matrix                                       
	double J_1_T[3][3];        	 // J^{-T}                                          
	double det_J;              // Jacobian (determinant of the Jacobian matrix)           
	double det_J_abs;          // Jacobian modulus                                        
	int type_of_hex;           // 0..30 - parallelepiped, 31..61 - hexahedron
	double grad_all[8][3];     // gradients from basis functions at integration point
	// local matrices 8x8 (regular, non-block) 
	double b[8][8];     // mass matrix    
	double hx, hy, hz;  // for the parallelepiped
	bool JforParCalc;   // whether the Jacobian was previously calculated in the case of a parallelepiped  
	void Calc_J(int n_of_point);
	void Calc_local_matrix_b_for_parallelepiped();
	void CalcMassMatrix();
	Hex_Local_Matrix(int i, long (*nver)[14], double (*xyz)[3]);
	~Hex_Local_Matrix();
};
