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
 *  This file contains headers for some basic routines for vector-matrix operations and error messages
 * 
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/COCR_FP/For_Solvers.h 
 *  Modified by D.Sc. Denis V. Vagin                                                                   
 *  Novosibirsk State Technical University,                                                            
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                   
 *  Version 2.0 December 10, 2024                                                                     
*/

#pragma once

double Scal(double *a, double *b, long n);    // Dot product 
double Norm_Euclid(double *a, long n);        // Euclidian norm of a vector
double Norm_Max(double *a, long n);           // Max norm of a vector 
double Projection_On_Axis(double *v,double *o); //Projection of the vector v onto the o axis
void Mult_Plot(double *a, double *x, double *y, long n); // Matrix-vector multiplication in dense format
void Mult_MV(long *ig, long *jg, double *ggl, double *ggu, double *di, double *x, double *y, long n); // Matrix-vector multiplication in sparse format 
double Relative_Error(double *analytic, double *numeric, long n); // Relative error
long Max_Long(long a, long b);   // Maximum of 2 numbers  
long Min_Long(long a, long b);   // Minimum of 2 numbers  
void Sort2(long *a, long *b);    // Sort 2 numbers        
double Interval(double *x, double *y);   // Distance between two 3D-points 
double Interval_Parallel_Lines(double *a0, double *a1, double *b0, double *b1);  // Distance between two parallel lines in 3D 
double Spline(double x, long n, double *xyz, double *values);                     // Linear interpolation
double Calc_dof(double *J, double *func, long n_local_edge);
void Memory_allocation_error(const char *var, const char *func);                 // reports a memory allocation error
void Cannot_open_file(const char *fname, const char *func);                      // reports a file open error and throws an exception 
void Cannot_open_file_but_continue(const char *fname, const char *func);         // reports an error opening the file, but the program continues
void Mult_Plot_AV(double *a, double *x, double *y, long n, long m);
