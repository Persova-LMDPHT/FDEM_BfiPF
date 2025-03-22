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
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/For_Solvers.h         
 *  Modified by D.Sc. Denis V. Vagin                                                                           
 *  Novosibirsk State Technical University,                                                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                          
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                           
 *  Version 2.0 December 10, 2024                                                                              
*/

#pragma once

double Scal(double *a, double *b, int n);
double Norm_Euclid(double *a, int n);
double Norm_Max(double *a, int n);
double Projection_On_Axis(double *v,double *o);
void Mult_Plot(double *a, double *x, double *y, int n);
void Mult_MV(int *ig, int *jg, double *ggl, double *ggu, double *di, double *x, double *y, int n);
double Relative_Error(double *analytic, double *numeric, int n);

void Memory_allocation_error(const char *var, const char *func);
void Cannot_open_file(const char *fname, const char *func);
void Cannot_open_file_but_continue(const char *fname, const char *func);
