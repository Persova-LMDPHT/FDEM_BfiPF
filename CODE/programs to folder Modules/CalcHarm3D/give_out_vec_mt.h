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
 *  This file contains headers for outputting E and EMF fields from solution in 3D VFEM
 *  
 *  
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/give_out_vec_mt.h     
 *  Modified by D.Sc. Denis V. Vagin                                                                           
 *  Novosibirsk State Technical University,                                                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                           
 *  Version 2.0 December 10, 2024                                                                              
*/

#pragma once
#include "AbstractFEM.h"
#include "OutputResultant3d.h"

class Give_out_vec_mt : public AbstractFEM3D
{
public:
	Vec_Prep_Data *d;
	T_Mapping_Vec *tmap;
	double *v3dat;
	OutputResultant3d *resultantA,*resultantB;
	std::complex <double> H_1d[2];
	std::complex <double> E_1d[2];
	std::complex <double> (*H)[3];
	std::complex <double> (*E)[3];
	double *impedance;
	double *rho;

	Give_out_vec_mt(Vec_Prep_Data *d, T_Mapping_Vec *tmap, double *v3dat);
	~Give_out_vec_mt();

	void Compute_1d_field();

	void Write_result_to_files();
	void Write_result_to_files_for_harm_loop();
	void Write_B_to_files_for_harm_loop(int StartType,vector<int> &RecvPlsIgB,int ParamInd);
	void Write_E_to_files_for_harm_loop(int StartType,vector<int> &RecvPlsIgE,int ParamInd);

	void Write_NULL_B_to_files_for_harm_loop(int ParamInd);
	void Write_NULL_E_to_files_for_harm_loop(int ParamInd);
	
	void Give_out_on_hex();

	int Read_1d_field(char *fname_in, double &Ex_s, double &Ex_c, double &Ey_s, double &Ey_c, 
		double &Hx_s, double &Hx_c,  double &Hy_s, double &Hy_c);

	int GetNumberOfNodes();
	int GetNumberOfElements();
	int GetElementNodesNumber();
	const pv::Point3D GetNode(const int& i_node);
	const pv::Point3D GetNodeTrue(const int& i_node);
	int GetNodeNumberOnElement(const int& i_element, const int& i_node);
	int GetElementMaterial(const int& i_element);
	int GetTypeOfElement(const int& i_element);
	double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type);
	int GetNumberOfResPoints(const Res3DValueType& r_type);
	pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point);
	int * GetPointerToRegular();
	int GetXSize();
	int GetYSize();
	int GetZSize();
	double *GetPointerToX();
	double *GetPointerToY();
	double *GetPointerToZ();
	void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& i_point, const int& i_time, int ipls);

	int ipls_cur;
	vector<Res3DValueType> vvta,vvtb;
};
