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
 *  This file contains headers for calculating nonstationary 3D VFEM task
 *  
 *  Written by D.Sc. Denis V. Vagin                                           
 *  Novosibirsk State Technical University,                                   
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                         
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)          
 *  Version 2.0 December 10, 2024                                              
 *  
*/

#pragma once
#include "T_Mapping.h"
#include "vec_prep_data.h"

class Time_approx_for_vfem 
{
public:

	int n;

	T_Mapping_Vec *Tmap;

	double finSP;

	int nmat;

	int *anomal_edges;
	int n_anomal_edges;

	Vec_Prep_Data *d;

	Time_approx_for_vfem(bool forPolygon, const void* _NormBFromPolygon, 
		bool flag_vfem_for_line, const void* _NormResFromLine);
	~Time_approx_for_vfem();

	int Read_data();

	// Finding and laying out centers and direction vectors of edges belonging to anomalous elements
	int Unload_anomal_nodes(char *fname);
	int Unload_anomal_nodesA0(char *fname);

	void GetRegularEdges(int i,vector<int> &re,int &size_re);

	vector<int> elemRenum;
	vector<bool> isElemAnomal;
	int nAnomalElem;

	int n_nodes_f;
	vector<int> isEdgesAnomal;
	vector<double> EdgesCrd[3];

	void AddWithRegularEdge(int iedge,vector<int> &isEdgesAnomal);

	int nzl;
	vector<double> zlay;
};
