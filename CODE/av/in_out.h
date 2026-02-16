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
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/in_out.h 
 *  Modified by D.Sc. Denis V. Vagin                                                               
 *  Novosibirsk State Technical University,                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                               
 * Version 2.0 December 10, 2024                                                                  
 */

#pragma once

class In_Out
{
private:
	int n_of_rec;
	int len_of_rec;
	int size_of_array;

	int type_of_values;

	char input_fname[100];
	char output_fname[100];

public:
	In_Out();
	~In_Out();

	int Write_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Write_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Read_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);
	int Read_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record);

	int Write_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record);
	int Write_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record);
	int Read_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record); 
	int Read_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record); 
	int Read_Bin_File_Of_Short(char *fname, short *massiv, int n_of_records, int len_of_record);

	int Read_int_From_Txt_File(char *fname, int *number);
	int Read_Double_From_Txt_File(char *fname, double *number);
	int Write_Double_To_Txt_File(char *fname, double number);
	int Write_int_To_Txt_File(char *fname, int number);

	int Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record);
	int Convert_File_Of_int_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record);

	int Read_inf2tr(char *fname, int *kuzlov, int *ktr, int *l1);
};
