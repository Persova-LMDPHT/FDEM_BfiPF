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
 *  This file contains code of basic routines for reading/writing files both in the binary and text formats
 *
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/in_out.cpp        
 *  Modified by D.Sc. Denis V. Vagin                                                                     
 *  Novosibirsk State Technical University,                                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                     
 *  Version 2.0 December 10, 2024                                                                         
*/

#include "stdafx.h"
#include "in_out.h"
#include "For_Solvers.h"
extern ofstream logfile;

In_Out::In_Out()
{
	len_of_rec = 0;
	n_of_rec = 0;
	size_of_array = 0;
}

In_Out::~In_Out()
{
}

int In_Out::Write_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record)
{
	FILE *fp;
	int i, j;
	
	if((fp=fopen(fname,"w"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_Txt_File_Of_int");
		return 1;
	}

	cout << "writing " << fname << "... ";
	logfile << "writing " << fname << "... ";
	for(i=0; i<n_of_records; i++)
	{
		for(j=0; j<len_of_record; j++)
			fprintf(fp,"%ld\t", massiv[i*len_of_record + j]+1);
		fprintf(fp,"\n");
	}
	logfile << " done\n";
	cout << " done\n";

	fclose(fp);
	return 0;
}

int In_Out::Write_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record)
{
	FILE *fp;
	int i, j;
	
	if((fp=fopen(fname,"w"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_Txt_File_Of_Double");
		return 1;
	}

	logfile << "writing " << fname;
	cout << "writing " << fname;

	fprintf(fp,"\n");
	fprintf(fp,"\n");

	for(i=0; i<n_of_records; i++)
	{
		for(j=0; j<len_of_record; j++)
			fprintf(fp,"%25.14e\t", massiv[i*len_of_record + j]);
		fprintf(fp,"\n");
	}

	logfile << " done\n";
	cout << " done\n";

	fclose(fp);
	return 0;
}

int In_Out::Write_kuslau_block(char *fname, int n_block, double eps, int maxiter)
{
FILE *fp;

	if((fp=fopen(fname,"w"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_kuslau_block");
		return 1;
	}
	fprintf(fp,"%ld\n", n_block);	
	fprintf(fp,"%ld\n", maxiter);	
	fprintf(fp,"%e\n", eps);
	fprintf(fp,"1e-16\n");

	fclose(fp);
	return 0;
}

int In_Out::Write_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record)
{
	int i, j;
	int temp;
	double value;
	FILE *fp;

	if((fp=fopen(fname,"w+b"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_Bin_File_Of_Double");
		return 1;
	}

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
		{
			value = massiv[i*len_of_record + j];
			temp = (int)fwrite(&value,sizeof(double),1,fp);
			if(temp!=1)
			{
				cout << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				logfile << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				fclose(fp);
				string str = "Cannot write to binary file ";
				str += fname;
				throw logic_error(str);
				return -1;
			}
		}
	fclose(fp);
	return 0;
}

int In_Out::Read_Txt_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record)
{
	FILE *fp;
	int i, j;
	
	if((fp=fopen(fname,"r"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Txt_File_Of_Double");
		return 1;
	}

	cout << "reading " << fname << "... ";
	logfile << "reading " << fname << "... ";

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
			fscanf(fp,"%lf", &massiv[i*len_of_record + j]);

	cout << " done\n";
	logfile << " done\n";

	fclose(fp);
	return 0;
}

int In_Out::Read_Txt_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record)
{
	FILE *fp;
	int i, j;
	
	if((fp=fopen(fname,"r"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Txt_File_Of_int");
		return 1;
	}

	logfile << "reading " << fname << "... ";
	cout << "reading " << fname << "... ";

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
			fscanf(fp,"%ld", &massiv[i*len_of_record + j]);

	logfile << " done\n";
	cout << " done\n";

	fclose(fp);
	return 0;
}

int In_Out::Read_Double_From_Txt_File(char *fname, double *number)
{
	FILE *fp;
	double temp;

	if((fp=fopen(fname,"r"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Double_From_Txt_File");
		return 1;
	}
	
	fscanf(fp,"%lf",&temp);
	*number = temp;

	fclose(fp);
	return 0;
}

int In_Out::Read_int_From_Txt_File(char *fname, int *number)
{
	FILE *fp;
	int temp;

	if((fp=fopen(fname,"r"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_int_From_Txt_File");
		return 1;
	}
	
	fscanf(fp,"%ld",&temp);
	*number = temp;

	fclose(fp);
	return 0;
}

int In_Out::Read_Bin_File_Of_Double(char *fname, double *massiv, int n_of_records, int len_of_record)
{
	int i, j;
	int temp;
	double value;
	FILE *fp;

	if((fp=fopen(fname,"r+b"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Bin_File_Of_Double");
		return 1;
	}

	logfile << "reading " << fname << "...\n";
	cout << "reading " << fname << "...\n";

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
		{
			temp = (int)fread(&value,sizeof(double),1,fp);
			if(temp!=1)	
			{
				cout << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				logfile << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				fclose(fp);

				string str = "Cannot read binary file ";
				str += fname;
				throw logic_error(str);

				return -1;
			}
			massiv[i*len_of_record + j] = value;
		}

	fclose(fp);
	return 0;
}

int In_Out::Read_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record)
{
	int i, j;
	int temp;
	int value;
	FILE *fp;

	if((fp=fopen(fname,"r+b"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Bin_File_Of_int");
		return 1;
	}

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
		{
			temp = (int)fread(&value,sizeof(int),1,fp);
			if(temp!=1)	
			{
				cout << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				logfile << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				fclose(fp);

				string str = "Cannot read binary file ";
				str += fname;
				throw logic_error(str);

				return -1;
			}
            massiv[i*len_of_record + j] = value;
		}

	fclose(fp);
	return 0;
}

int In_Out::Read_Bin_File_Of_Short(char *fname, short *massiv, int n_of_records, int len_of_record)
{
	int i, j;
	int temp;
	short value;
	FILE *fp;

	if((fp=fopen(fname,"r+b"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Read_Bin_File_Of_int");
		return 1;
	}

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
		{
			temp = (int)fread(&value,sizeof(short),1,fp);
			if(temp!=1)	
			{
				cout << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				logfile << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				fclose(fp);

				string str = "Cannot read binary file ";
				str += fname;
				throw logic_error(str);

				return -1;
			}
            massiv[i*len_of_record + j] = value;
		}

	fclose(fp);
	return 0;
}

int In_Out::Write_Bin_File_Of_Long(char *fname, int *massiv, int n_of_records, int len_of_record)
{
	int i, j;
	int temp;
	int value;
	FILE *fp;

	if((fp=fopen(fname,"w+b"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_Bin_File_Of_int");
		return 1;
	}

	for(i=0; i<n_of_records; i++)
		for(j=0; j<len_of_record; j++)
		{
			value = massiv[i*len_of_record + j]+1;
			temp = (int)fwrite(&value,sizeof(int),1,fp);
			if(temp!=1)
			{
				cout << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';
				logfile << "feof=" << feof(fp) << " ferror=" << ferror(fp) << '\n';

				string str = "Cannot write to binary file ";
				str += fname;
				throw logic_error(str);

				fclose(fp);
				return -1;
			}
        }

	fclose(fp);
	return 0;
}

int In_Out::Convert_File_Of_Double_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record)
{
	double *array_of_double=NULL;

	size_of_array = n_of_records*len_of_record;

	array_of_double = new double[size_of_array];
	if(array_of_double == 0)
		Memory_allocation_error("array_of_double", "Convert_File_Of_Double_From_Bin_To_Txt");

	Read_Bin_File_Of_Double(file_in, array_of_double, n_of_records, len_of_record);

	Write_Txt_File_Of_Double(file_out, array_of_double, n_of_records, len_of_record);

	if(array_of_double) {delete [] array_of_double; array_of_double=NULL;}

	return 0;
}

int In_Out::Convert_File_Of_int_From_Bin_To_Txt(char *file_in, char *file_out, int n_of_records, int len_of_record)
{
	int *array_of_int=NULL;

	size_of_array = n_of_records*len_of_record;

	array_of_int = new int[size_of_array];
	if(array_of_int == 0)
		Memory_allocation_error("array_of_int", "In_Out::Convert_File_Of_int_From_Bin_To_Txt");

	Read_Bin_File_Of_Long(file_in, array_of_int, n_of_records, len_of_record);

	Write_Txt_File_Of_Long(file_out, array_of_int, n_of_records, len_of_record);

	if(array_of_int) {delete [] array_of_int; array_of_int=NULL;}

	return 0;
}

int In_Out::Write_jg(char *fname, int *ig, int *jg, int n)
{
	FILE *fp;
	int i, j;
	
	if((fp=fopen(fname,"w"))==0)
	{
		Cannot_open_file_but_continue(fname, "In_Out::Write_jg");
		return 1;
	}

	cout << "writing jg...";
	logfile << "writing jg...";
	for(i=0; i<n; i++)
	{
		for(j=ig[i]; j<ig[i+1]; j++)
			fprintf(fp,"%ld\t", jg[j] + 1);
		fprintf(fp,"\n");
	}
	printf("done\n");
	logfile << " done\n";

	fclose(fp);
	return 0;
}

int In_Out::Write_kuslau(char *fname, int n, double eps, int maxiter)
{
	FILE *fp;

	if((fp=fopen(fname,"w"))==0)
	{
		printf("Error: Cannot open file \"%s\" for writing.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for writing.\n";
		return 1;
	}
	fprintf(fp,"%ld\n", n);	
	fprintf(fp,"%e\n", eps);	
	fprintf(fp,"%ld\n", maxiter);	

	fclose(fp);
	return 0;
}

int In_Out::Menu(char *input_fname, char *output_fname)
{
	int c;
	int len_of_rec;
	int n_of_rec;

	strcpy(this->input_fname, input_fname);
	strcpy(this->output_fname, output_fname);

	printf("Enter type of variables:\n");
	printf("1 - int\n");
	printf("2 - double\n");
	scanf("%ld", &c);

	printf("Enter number of records:\n");
	scanf("%ld", &n_of_rec);

	printf("Enter lenght of record:\n");
	scanf("%ld", &len_of_rec);

	switch(c)
	{
	case 1:
		this->Convert_File_Of_int_From_Bin_To_Txt(this->input_fname,this->output_fname,n_of_rec,len_of_rec);
		break;
	case 2:
		this->Convert_File_Of_Double_From_Bin_To_Txt(this->input_fname,this->output_fname,n_of_rec,len_of_rec);
	}

	return 0;
}

int In_Out::Read_inftry(char *fname, int *kuzlov, int *kpar, int *l1)
{
	char buffer[100];
	FILE *fp;
	
	if((fp=fopen(fname,"r"))==0)
	{
		printf("Error: Cannot open file \"%s\" for reading.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for reading.\n";
		return 1;
	}

	while(!feof(fp))
	{
		fscanf(fp, "%s", buffer);
		if(strcmp(buffer,"KUZLOV=")==0)
		{
			fscanf(fp, "%ld", kuzlov);
			continue;
		}
		if(strcmp(buffer,"KPAR=")==0)
		{
			fscanf(fp, "%ld", kpar);
			continue;
		}
		if(strcmp(buffer,"KT1=")==0)
		{
			fscanf(fp, "%ld", l1);
			continue;
		}
	}

	fclose(fp);
	return 0;
}

int In_Out::Write_Double_To_Txt_File(char *fname, double number)
{
	FILE *fp;

	if((fp=fopen(fname,"w"))==0)
	{
		printf("Error: Cannot open file \"%s\" for writing.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for writing.\n";
		return 1;
	}
	
	fprintf(fp,"%25.13e",number);

	fclose(fp);
	return 0;
}

int In_Out::Write_int_To_Txt_File(char *fname, int number)
{
	FILE *fp;

	if((fp=fopen(fname,"w"))==0)
	{
		printf("Error: Cannot open file \"%s\" for writing.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for writing.\n";
		return 1;
	}

	fprintf(fp,"%ld",number);

	fclose(fp);
	return 0;
}

int In_Out::Read_n_eps_maxiter(char *fname, int *n, double *eps, int *maxiter)
{
	FILE *fp;
	int t1;
	double t2;
	int t3;

	if((fp=fopen(fname, "r"))==0)
	{
		printf("Error: Cannot open file \"%s\" for reading.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for reading.\n";
		return 1;
	}
	fscanf(fp, "%ld", &t1);
	fscanf(fp, "%lf", &t2);
	fscanf(fp, "%ld", &t3);

	*n = t1;
	*eps = t2;
	*maxiter = t3;

	fclose(fp);
	return 0;
}

int In_Out::Read_n_maxiter_eps(char *fname, int *n, int *maxiter, double *eps)
{
	FILE *fp;

	if((fp=fopen(fname, "r"))==0)
	{
		printf("Error: Cannot open file \"%s\" for reading.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for reading.\n";
		return 1;
	}
	fscanf(fp, "%ld", n);
	fscanf(fp, "%ld", maxiter);
	fscanf(fp, "%lf", eps);	

	fclose(fp);
	return 0;
}

int In_Out::Read_inf2tr(char *fname, int *kuzlov, int *ktr, int *l1)
{
	char buffer[100];
	FILE *fp;

	if((fp=fopen(fname,"r"))==0)
	{
		printf("Error: Cannot open file \"%s\" for reading.\n",fname);
		logfile << "Error: Cannot open file " << fname << " for reading.\n";
		return 1;
	}

		while(!feof(fp))
	{
		fscanf(fp, "%s", buffer);
		if(strcmp(buffer,"KUZLOV=")==0 || strcmp(buffer,"kuzlov=")==0)
		{
			fscanf(fp, "%ld", kuzlov);
			continue;
		}
		if(strcmp(buffer,"KTR=")==0 || strcmp(buffer,"ktr=")==0)
		{
			fscanf(fp, "%ld", ktr);
			continue;
		}
		if(strcmp(buffer,"KT1=")==0 || strcmp(buffer,"kt1=")==0)
		{
			fscanf(fp, "%ld", l1);
			continue;
		}
	}

	fclose(fp);
	return 0;
}

