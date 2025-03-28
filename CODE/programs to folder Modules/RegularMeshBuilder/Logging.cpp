/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *                    
 * The file contains basic utility functions and simple operations
 *
 *  Earlier uploaded: https://github.com/kiselev2013/A-E_Formulations_3D_TimeDomain_EM_HorizontalGroundedWireSource/blob/master/CODE/programs%20to%20folder%20Modules/RegularMeshBuilder/Logging.cpp
 *  Novosibirsk State Technical University,                                                                                                                                                       
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                                                             
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                                                              
 *  Version 2.0 16 January, 2023                                                                                                                                                                  
*/

#include "Logging.h"

FILE *log_file;

// open file for writing
int open_file_w(char *file_name, FILE **file_stream)
{
	char buf[2048];

	if (!((*file_stream) = fopen(file_name, "w")))
	{
		sprintf(buf, "Error : Could not write file '%s'\n", file_name);
		write_to_log(buf);
		return 1;
	}

	return 0;
}

// open log file
int open_log(char *file_name)
{
	if (open_file_w(file_name, &log_file) != 0)
	{
		printf("Error : could not open file '%s'\n", file_name);
		return 1;
	}

	return 0;
}

// log message
void write_to_log(char *str)
{
	printf("%s", str);
	fprintf(log_file, "%s", str);
	fflush(log_file);
}

// log message
void write_to_log(const char *str)
{
	printf("%s", str);
	fprintf(log_file, "%s", str);
	fflush(log_file);
}