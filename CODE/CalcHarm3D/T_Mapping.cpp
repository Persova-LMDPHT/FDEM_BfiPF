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
 *  This file contains code for working with edge mesh in 3D VFEM
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 * Version 2.0 December 10, 2024
 *  
*/

#include "stdafx.h"
#include "T_Mapping.h"
#include "Portret.h"
extern ofstream logfile;

T_Mapping_Vec::T_Mapping_Vec(long (*nver)[14], double (*xyz)[3], long kuzlov, long kpar)
{
	this->nver = nver;
	this->xyz = xyz;
	this->kuzlov = kuzlov;
	this->kpar = kpar;

	ed = NULL;
	edges = NULL;
}

T_Mapping_Vec::~T_Mapping_Vec()
{
	if(ed) {delete [] ed; ed = NULL;}
	if(edges) {delete [] edges; edges = NULL;}
}
