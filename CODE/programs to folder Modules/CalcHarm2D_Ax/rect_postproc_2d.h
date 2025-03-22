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
 *  This file contains headers for structures for finite element mesh storing
 * 
 *  Written by D.Sc. Denis V. Vagin                                         
 *  Novosibirsk State Technical University,                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 December 10, 2024                                           
*/

#pragma once

class Rect_preproc_data
{
public:
	Rect_preproc_data();
	~Rect_preproc_data();
	int Read_rect_preproc_data();
	int Read_mesh_rz();
	int Read_xyzVectorE0();
	int Read_xyzVectorB();
	int Read_xyzVectorE();
	int Read_xyzVectorB2d();
	int Read_xyzVectorE2d();
	int Read_pointres();
	int Read_1d();
	long Rectangle_or_quadrilateral(long num, double (*xy)[2], long (*nver)[6]);
	enum PrTypeForLine {ptNone, ptAphi, ptHphi, ptAr, ptAphip} pType;
	int alfa;
	long n_materials;
	double *mu2d;
	double *mu0;
	double nu;
	double *sigma2d;
	double *sigma0;
	double *dpr2d;
	long tsize2d;
	long ig_t_n_1;
	long kuzlov;
	long n_nodes_c;
	long krect;
	long kt1;
	long *l1;
	short *nvk1;
	long (*nvtr)[4];
	long *type_of_rect;
	long *nvkat;
	double (*xy)[2];
	double *pointres;
	long n_pointres;
	long n_1d;
	double *coords_1d;
	double *sin_1d;
	double *cos_1d;
	int *mtr3d2d;
	long n_xyzVectorE0;
	double (*xyzVectorE0)[3];
	long n_xyzVectorB;
	double (*xyzVectorB)[3];
	long n_xyzVectorE;
	double (*xyzVectorE)[3];
	long n_xyzVectorB2d;
	double (*xyzVectorB2d)[3];
	long n_xyzVectorE2d;
	double (*xyzVectorE2d)[3];
	vector<long> reg;
	long qr;
	long qz;
	vector<double> rm;
	vector<double> zm;
	long nreg;
};

