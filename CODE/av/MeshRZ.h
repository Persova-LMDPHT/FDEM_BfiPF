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
 *  Written by  D.Sc. Denis V. Vagin                                        
 *  Novosibirsk State Technical University,                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 December 10, 2024                                           
*/

#pragma once

class MeshRZ
{
public:

	int kuzlov;
	int n_nodes_c;
	int krect;
	int kt1;
	int *l1;

	int n_materials;
	double *sigma;
	vector<double> sigma_xy;
	vector<double> sigma_z;
	double *mu;

	int tsize2d;
	double (*xy)[2];
	int (*nvtr)[4];
	int *type_of_rect;
	int *nvkat;

	int nreg, qr, qz;
	int *reg;
	double *rm, *zm;

	int ntime;
	double *time;

	int n_2kr;
	int (*edges_2kr)[2];
	vector <int> nodes_2kr;

	double az_line;
	double bz_line;

	int ReadMeshRZ(char *path);
	int Rectangle_or_quadrilateral(int num, double (*xy)[2], int (*nver)[6]);
	int Read2Kr(char *f_n, char *f_2kr);
	int Find2kr(double az, double bz, double r0, double eps);
	int Find1kr(double eps);
	int Read_Az_Bz(char *fname);

	int FindNodeWithDeltaFunction(double z0, int *elem1, int *elem2);
	void FixBoundCond();

	MeshRZ();
	~MeshRZ();
};
