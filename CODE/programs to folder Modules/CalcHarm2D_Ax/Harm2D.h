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
 *  This file contains headers of function for calculating the eleptic source part of field for harmonic task
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                       
 *  Novosibirsk State Technical University,                                                          
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)      
 * Version 2.0 December 10, 2024                                                                                                                       
*/

#pragma once
#include "rect_postproc_2d.h"

int SolveHarmonic2D(Rect_preproc_data *d,double *v2, double cd, int npls,int nfreq);
int ProcTask2DLine_Harmonic(double curr, int nfreq);
int FindNodeWithDeltaFunction(long kuzlov, double (*rz)[2], double rd, double zd, long *node);

