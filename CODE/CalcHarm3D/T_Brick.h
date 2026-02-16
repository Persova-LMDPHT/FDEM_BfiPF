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
 *  This file contains structures and headers for some basic routines for generation of local matrices and local vector of the right-hand side.
 *  Calculation of values of basis functions inside a finite element.
 *  The T_Brick class contains functions for working with vector (edge-elements) for hexahedron mesh:
 * 	- calculation of local stiffness and mass matrices;
 * 	- calculation of the vector of the right-hand side through E_normal;
 *	- calculation of the Jacobi matrix at Gauss points and at an arbitrary point inside the hexahedron;
 * 	- output a solution and a curl inside a hexahedron or parallelepiped;
 * Nodal basis functions are needed here to calculate the transformation.
 * Template element [-1, 1]^3. (for both nodal and vector).
 * Vector basis functions with tangential components 2/hx, 2/hy, 2/hz along edges.
 *
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/T_Brick.h    
 *  Modified by D.Sc. Denis V. Vagin                                                                  
 *  Novosibirsk State Technical University,                                                           
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                 
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                  
 *  Version 2.0 December 10, 2024                                                                      
*/

#pragma once
#include "vec_prep_data.h"

class T_Brick
{
public:
	double hx, hy, hz;
	double xk, xk1, yk, yk1, zk, zk1;
    long num;

	double b[12][12];
	double c[12][12];
	double c0[12][12];
	double cSigma[12][12];
	double cSigma0[12][12];
	double cDpr[12][12];
	double cDpr0[12][12];

	double f_re[12];

	double a[12][12];
	double g[12];
	double g_harm[24];

	double g_re_sig[12], g_im_sig[12];
	double g_re_dpr[12], g_im_dpr[12];
	double g_re_b[12], g_im_b[12];
	double g8[8];

	double mu;
	double mu0;
	double sigma;
    Tensor sigmaTensor, sigmaTensor0, s_s0, dprTensor, dprTensor0, d_d0;
	double sigma0;
	double dpr;
	double dpr0;
	long n_mat;

	long (*nver)[14];
	long (*ed)[25];
	long *nvkat;
	long (*edges)[2];
	double (*xyz)[3];

	double *En;
	double (*En_nodes)[3];

	double mrot[3][3];

	long (*nvetr)[20];
	double MtrMass[20][20];
	double MtrMassRight[20][12];
	double MtrGest[20][20];
	double VctRight[40];
	double g_re_sig_av[20], g_im_sig_av[20];

	int tasktype;

	T_Brick(double *x_coords, double *y_coords, double *z_coords, long type_of_hex);

	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double *En);

	T_Brick(long num, long (*nver)[14], long (*ed)[25], long (*edges)[2], double (*xyz)[3], long *nvkat,
		double *sigma3d, double *sigma0, double *mu3d, double omega,
		long alpha, long n_1d, double *z_1d, double *sin_1d, double *cos_1d);

	T_Brick(long num, long (*nver)[14], double (*xyz)[3]);

	T_Brick();

	~T_Brick();

	void Compute_Local_Matrix_And_Vector(const long what_compute); 
	// Calculate the local stiffness matrix of the element
	void Compute_Local_Matrix_B(); 
	// Calculate the local mass matrix of an element
	void Compute_Local_Matrix_C(); 
	void Compute_Local_Matrix_C_Tensor(Tensor &t,double mm[][12]);
	void Compute_Local_Vector_For_Anomal_Problem();
	void ComputeLocalVectorMuEpsSigma(double *An, double *d2An);
	void GetRotationMatrix();
	void RotateTensors();

	double omega;
	double asin0[12], acos0[12];
	long alpha;
	long n_1d;
	double *z_1d;
	double *sin_1d;
	double *cos_1d;

	void Calc_local_vector_for_MT();
	void Calc_asin_acos_at_middle_of_edges();

	void Calc_block_local_matrix_and_vector(bool fclcmtrx,bool fclcrprt);
	void Calc_J_Node_2(double x, double y, double z);
	double dPhi_node_2(long i, long j, double x, double y, double z);
	void Calc_V_Node_2(double *q,double x, double y, double z);

	double V[3];
	double x[8], y[8], z[8];
	double J[3][3];
	double J_1[3][3];
	double J_1_T[3][3];
	double det_J;
	double det_J_abs;
	long type_of_hex;
	double phi_all[12][3];
	double rot_all[12][3];

	void Mapping(double *in, double *out);

	void Calc_J(int n_of_point);
	void Calc_J(double x, double y, double z);
	void Calc_J_in_parallelepiped();

	void Calc_local_matrix_b_for_hexahedron();
	void Calc_local_matrix_c_for_hexahedron();

	void Calc_value_inside_hex(double *ves, double *in, double *out); 
	void Basis_func_on_vec_par(long i, double ves, double *in, double *out);
	void Basis_func_on_reference_vec_par(long i, double *in, double *out);	
	void Basis_func_on_vec_hex(long i, double ves, double *in, double *out);
	void Calc_rotor_inside_hex(double *ves, double *in, double *out); 
	void Calc_rotor_inside_hex_crd(double *ves, double *in, double *out,int ind,
		double *cfbx, double *cfby, double *cfbz); 
	void Calc_rotor_coefficients_inside_hex_crd(double *in, int ind, double *cfbx, double *cfby, double *cfbz);
	void Rot_of_basis_func_on_reference_vec_par(long i, double *in, double *out);
	void Rotx_of_basis_func_on_reference_vec_par(long i, double x, double *out);
	void Roty_of_basis_func_on_reference_vec_par(long i, double y, double *out);
	void Rotz_of_basis_func_on_reference_vec_par(long i, double z, double *out);
	void Rot_of_basis_func_on_vec_par(long i, double ves, double *in, double *out);
	void Rot_of_basis_func_on_vec_par_x(int i, double ves, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_y(int i, double ves, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_z(int i, double ves, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_x(int i, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_y(int i, double *in, double &out);
	void Rot_of_basis_func_on_vec_par_z(int i, double *in, double &out);

	double l0(double x);
	double l1(double x);

	double Phi_node(long i, double x, double y, double z);

	double dPhi_node(long i, long j, double x, double y, double z);

	void VectorFieldOnPar(double x, double y, double z, double *ves,
		double *x_out, double *y_out, double *z_out);
	double ScalarFieldOnPar(double x, double y, double z, double *ves);
	double DxOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DyOfScalarFieldOnPar(double x, double y, double z, double *ves);
	double DzOfScalarFieldOnPar(double x, double y, double z, double *ves);

	void ScalarFieldOnParCff(double x, double y, double z, double *cff);

	double GetValueInHexCenter(double *q);
	void GetGradInHexCenter(double *q, double *out, int cj);
	void VectorFieldXOnPar3(double y, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);
	void VectorFieldYOnPar3(double x, double z, double *ves_j2, double *ves_j1, double *ves_j,
		double *out_j2, double *out_j1, double *out_j);

	void RotXOnPar(double x, double *ves, double *out, bool loc_c=false);
	void RotYOnPar(double y, double *ves, double *out, bool loc_c=false);
	void RotZOnPar(double z, double *ves, double *out);
	void RotZOnPar3(double z,
		double *ves1, double *ves2, double *ves3,
		double *out1, double *out2, double *out3);

	void Transformation_of_variables(double *in, double *out);
	void Transformation_of_variables(double *x, double *y, double *z);
	double Xi(double x);
	double Eta(double y);
	double Zeta(double z);

	void Set_dpr(double dpr);
	void Set_dpr0(double dpr0);
	void Set_mu0(double mu0);
	
	Vec_Prep_Data *d;


	double asin0n[8][3], acos0n[8][3];
	double asin0c[3], acos0c[3];

	void Calc_asin_acos_at_nodes();

	void GetVectorFieldNodes(double *ves, double *ax, double *ay, double *az);

	int npls,ipls,n_edges;

	void Calc_ss0();

	void CopyMaterial(int n_mat, double *sigma3d, double *sigma0, double *mu3d);
};

const double MIDDLE_OF_LOCAL_EDGE[12][3] = {
	 0.0, -1.0, -1.0,
	 0.0,  1.0, -1.0,
	 0.0, -1.0,  1.0,
	 0.0,  1.0,  1.0,

	-1.0,  0.0, -1.0,
	-1.0,  0.0,  1.0,
	 1.0,  0.0, -1.0,
	 1.0,  0.0,  1.0,

	-1.0, -1.0,  0.0,
	 1.0, -1.0,  0.0,
	-1.0,  1.0,  0.0,
	 1.0,  1.0,  0.0
};

const double LOCAL_COORDS_OF_NODES[8][3] = 
{
	-1.0, -1.0, -1.0,
	 1.0, -1.0, -1.0,
	-1.0,  1.0, -1.0,
	 1.0,  1.0, -1.0,

	-1.0, -1.0,  1.0,
	 1.0, -1.0,  1.0,
	-1.0,  1.0,  1.0,
	 1.0,  1.0,  1.0
};

const double TANGENT_VECTORS_ON_REFERENCE_CUBE[12][3] = {
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,
	1.0, 0.0, 0.0,

	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 1.0, 0.0,

	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0,
	0.0, 0.0, 1.0
};

const long REG_EDGES[12][2]={0,1, 2,3, 4,5, 6,7, 0,2, 4,6, 1,3, 5,7, 0,4, 1,5, 2,6, 3,7 };
