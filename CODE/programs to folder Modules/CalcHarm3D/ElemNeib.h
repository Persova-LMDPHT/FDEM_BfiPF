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
 *  This file contains headers of functions for reading and writing a structure of elements neighbors
 *
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024  
*/

#pragma once

struct ElemNeib{
	vector<int> neib[6];  // lists of neighbors for element surfaces
	ElemNeib(){
		neib[0].clear();neib[1].clear();neib[2].clear();
		neib[3].clear();neib[4].clear();neib[5].clear();
	}
};	

int ReadElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);
int WriteElemNeib(vector<ElemNeib> &ElemNeibVec,int kpar);
