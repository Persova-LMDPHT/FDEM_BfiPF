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
 *  This file contains headers for direct FEM output E and B functions.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#pragma once
#include "task2d.h"

struct RectOutData
{
	double sigma;
	double hr, hz, psi, eta, psi2, eta2;
	int node0, node1, node2, node3;
};

void Output_Field(int nrec, double *v2es, double *v2ec, double *results, double *resultc, RectOutData *vRecOutData);
int FindElementsForReceivers(int nrec, int qr, int qz, double *rm, double *zm, int *reg, PointRZ *rec, RectOutData *vRecOutData,
	PointRZ *pnt, Rect *rect, double *sigma, double rmin, double zmin, double rmax, double zmax);
