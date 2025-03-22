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
 *  This file contains headers for converting matrix from block CSRC to CSR format
 *  
 *  
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/FormatConverter.h  
 *  Novosibirsk State Technical University,                                                                  
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                         
 *  Version 1.2 April 7, 2021                                                                                
 *  
*/

#pragma once

class FormatConverter
{
public:
	void FromRSFToCSR_Real_1_Sym(int nb, int *ig, int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_Sym(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
	void From2x2ToCSR_Complex_1_NS(int nb, int *ig, int *idi, int *ijg,int *sz_iptr, int *sz_jptr);
		void From2x2ToCSRComplex_2_Sym(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block,
			MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void FromRSFToCSR_Real_2_Sym(int nb, int *ig, int *jg, double *di, double *gg,MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	void From2x2ToCSRComplex_2_NS(int nb, int *ig, int *jg, int *idi, int *ijg,double *di_block, double *ggl_block, double *ggu_block,
		MKL_INT64 *iptr, MKL_INT64 *jptr, double *aelem);
	FormatConverter();
	~FormatConverter();
};
