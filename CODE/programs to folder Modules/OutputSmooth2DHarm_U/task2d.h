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
 *  This file contains headers for mesh structures and functions
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#pragma once

struct PointRZ
{
	double r, z;
};

struct Rect
{
	int nodes[5];
	int rtype;
	int mtr;
};

struct task2d
{
	char path[256];

	int kpnt, krect, nc;
	PointRZ* pnt;
	Rect* rect;

	int nmat3d;
	int *mtr3d2d;

	int nmat2d;
	double *sigma,*sigmaZ;
	double *dpr;

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;

	task2d(char *_path);
	~task2d();

	int Read();
};
