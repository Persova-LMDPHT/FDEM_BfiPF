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
 *  This file contains headers for basic routines for reading/writing files both in the binary and text formats
 *
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/in_out.h          
 *  Novosibirsk State Technical University,                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                      
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                       
 *
*/ 

#pragma once

class In_Out
{
private:
	long n_of_rec;
	long len_of_rec;
	long size_of_array;
	long type_of_values;
	char input_fname[100];
	char output_fname[100];
public:
	In_Out();
	~In_Out();
	int Write_Txt_File_Of_Double(char *fname, double *massiv, long n_of_records, long len_of_record);
	int Write_Bin_File_Of_Double(char *fname, double *massiv, long n_of_records, long len_of_record);
	int Read_Txt_File_Of_Double(char *fname, double *massiv, long n_of_records, long len_of_record);
	int Read_Bin_File_Of_Double(char *fname, double *massiv, long n_of_records, long len_of_record);
	int Write_Txt_File_Of_Long(char *fname, long *massiv, long n_of_records, long len_of_record);
	int Write_Bin_File_Of_Long(char *fname, long *massiv, long n_of_records, long len_of_record);
	int Read_Txt_File_Of_Long(char *fname, long *massiv, long n_of_records, long len_of_record); 
	int Read_Bin_File_Of_Long(char *fname, long *massiv, long n_of_records, long len_of_record); 
	int Read_Bin_File_Of_Short(char *fname, short *massiv, long n_of_records, long len_of_record);
	int Read_Long_From_Txt_File(char *fname, long *number);
	int Read_Double_From_Txt_File(char *fname, double *number);
	int Write_Double_To_Txt_File(char *fname, double number);
	int Write_Long_To_Txt_File(char *fname, long number);
	int Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out, long n_of_records, long len_of_record);
	int Convert_File_Of_Long_From_Bin_To_Txt(char *file_in, char *file_out, long n_of_records, long len_of_record);
	int Menu(char *input_fname, char *output_fname);
	int Write_jg(char *fname, long *ig, long *jg, long n);
    int Write_kuslau(char *fname, long n, double eps, long maxiter);
	int Write_kuslau_block(char *fname, long n_block, double eps, long maxiter);
	int Read_inftry(char *fname, long *kuzlov, long *kpar, long *l1);
	int Write_inftry(char *fname, long kuzlov, long kpar, long l1);
	int Read_1d_data(long n_1d, double *coords_1d, double *sin_1d, double *cos_1d);
	int Read_n_eps_maxiter(char *fname, long *n, double *eps, long *maxiter);
	int Read_n_maxiter_eps(char *fname, long *n, long *maxiter, double *eps);
	int Read_inf2tr(char *fname, long *kuzlov, long *ktr, long *l1);
	int Write_inf2tr(char *fname, long kuzlov, long ktr, long l1);
};
