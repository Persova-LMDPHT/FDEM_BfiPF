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
 *  This file contains main function for calculating the grounded source part of field for harmonic task
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                             
 *  Novosibirsk State Technical University,                                                                
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                      
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)            
 * Version 2.0 December 10, 2024                                                
*/

#include "stdafx.h"
#include "VelHarm2d.h"
#include "MeshRZ.h"
#include "Portret.h"
#include "global_slae_2d.h"
#include "in_out.h"
#include "bound_cond_2d.h"

ofstream logfile;

int _tmain(int argc, TCHAR* argv[], TCHAR* envp[])
{
	int nRetCode = 0;
	{
		try
		{
			OutputB();			
		}
		catch(...)
		{
			nRetCode = 1;
		}
	}

	return nRetCode;
}
