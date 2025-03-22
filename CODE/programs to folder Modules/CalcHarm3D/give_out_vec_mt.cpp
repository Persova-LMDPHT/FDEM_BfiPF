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
 *  This file contains code for outputting E and B fields from solution in 3D VFEM
 *  
 *  Based version: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/give_out_vec_mt.cpp
 *  Modified by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 *  Version 2.0 December 10, 2024
 *  
*/

#include "stdafx.h"
#include "OutputArbitrary.h"
#include "give_out_vec_mt.h"
extern ofstream logfile;
Give_out_vec_mt::Give_out_vec_mt(Vec_Prep_Data *d, T_Mapping_Vec *tmap, double *v3dat)
{
	this->d = d;
	this->tmap = tmap;
	this->v3dat = v3dat;

	this->H = new std::complex<double>[d->n_pointresB*Subdomain::npls][3];
	this->E = new std::complex<double>[d->n_pointresE*Subdomain::npls][3];
	this->impedance = new double[d->n_pointresB*Subdomain::npls];
	this->rho = new double[d->n_pointresB*Subdomain::npls];

	if (d->tasktype==0)
		Compute_1d_field();
}
Give_out_vec_mt::~Give_out_vec_mt()
{
	delete [] H;
	delete [] E;
	delete [] impedance;
	delete [] rho;
}
void Give_out_vec_mt::Compute_1d_field()
{
	double usin, ucos, ducos, dusin;
	double alpha, w, mu;

	w = d->nu*2.0*PI;
	alpha = d->alfa;
	mu = MU_0;
	
	usin = d->usin[d->n_1d-1];
	ucos = d->ucos[d->n_1d-1];

	ducos = 0.0;
	dusin = mu;

	E_1d[0] = std::complex<double>(w*alpha*ucos, -w*alpha*usin);
	E_1d[1] = std::complex<double>(w*(1.0 - alpha)*ucos, w*(alpha - 1.0)*usin);

	H_1d[0] = std::complex<double>((alpha - 1.0)*dusin/mu, (alpha - 1.0)*ducos/mu);
	H_1d[1] = std::complex<double>(alpha*dusin/mu, alpha*ducos/mu);
}
void Give_out_vec_mt::Give_out_on_hex()
{
	int i, j;
	double w;
	const double mu = MU_0;
	w = d->nu*2.0*PI;
	if(d->n_pointresB)
	{
		resultantB->ValueType=vtRotzASin;
		resultantB->Output(0,ipls_cur,0,vvtb);
		resultantB->ValueType=vtRotzACos;
		resultantB->Output(0,ipls_cur,1,vvtb);
		resultantB->ValueType=vtRotxASin;
		resultantB->Output(0,ipls_cur,2,vvtb);
		resultantB->ValueType=vtRotxACos;
		resultantB->Output(0,ipls_cur,3,vvtb);
		resultantB->ValueType=vtRotyASin;
		resultantB->Output(0,ipls_cur,4,vvtb);
		resultantB->ValueType=vtRotyACos;
		resultantB->Output(0,ipls_cur,5,vvtb);
	}
	if(d->n_pointresE)
	{
		resultantA->ValueType=vtAzSin;
		resultantA->Output(0,ipls_cur,0,vvta);
		resultantA->ValueType=vtAzCos;
		resultantA->Output(0,ipls_cur,1,vvta);
		resultantA->ValueType=vtAxSin;
		resultantA->Output(0,ipls_cur,2,vvta);
		resultantA->ValueType=vtAxCos;
		resultantA->Output(0,ipls_cur,3,vvta);
		resultantA->ValueType=vtAySin;
		resultantA->Output(0,ipls_cur,4,vvta);
		resultantA->ValueType=vtAyCos;
		resultantA->Output(0,ipls_cur,5,vvta);
	}
}
void Give_out_vec_mt::Write_result_to_files()
{
	long j;
	FILE *ex_c, *ex_s, *ey_c, *ey_s; 
	FILE *hx_c, *hx_s, *hy_c, *hy_s, *hz_c, *hz_s; 
	FILE *rok;
	FILE *f_impedance;

	ofstream out;

	out.open("normal_field");

	out << "Ex_s:\t" << real(E_1d[0]) << endl;
	out << "Ex_c:\t" << imag(E_1d[0]) << endl;
	out << "Ey_s:\t" << real(E_1d[1]) << endl;
	out << "Ey_c:\t" << imag(E_1d[1]) << endl;

	out << "Hx_s:\t" << real(H_1d[0]) << endl;
	out << "Hx_c:\t" << imag(H_1d[0]) << endl;
	out << "Hy_s:\t" << real(H_1d[1]) << endl;
	out << "Hy_c:\t" << imag(H_1d[1]) << endl;

    out.close();
	out.clear();

	ex_c = fopen("ex_c", "w"); ex_s = fopen("ex_s", "w");
	ey_c = fopen("ey_c", "w"); ey_s = fopen("ey_s", "w");
	hx_c = fopen("hx_c", "w"); hx_s = fopen("hx_s", "w");
	hy_c = fopen("hy_c", "w"); hy_s = fopen("hy_s", "w");
	hz_c = fopen("hz_c", "w"); hz_s = fopen("hz_s", "w");
	rok = fopen("rok", "w");
	f_impedance = fopen("impedance", "w");

	fprintf(ex_c, "\n\n");
	fprintf(ex_s, "\n\n");
	fprintf(ey_c, "\n\n");
	fprintf(ey_s, "\n\n");
	fprintf(hx_c, "\n\n");
	fprintf(hx_s, "\n\n");
	fprintf(hy_c, "\n\n");
	fprintf(hy_s, "\n\n");	
	fprintf(hz_c, "\n\n");
	fprintf(hz_s, "\n\n");
	fprintf(rok, "\n\n");
	fprintf(f_impedance, "\n\n");

	for(j=0; j<d->n_pointresB; j++)
	{
		double x;
		x = d->pointresB[j][0];
		fprintf(ex_c, "%g %g\n", x, imag(E[j][0]));
		fprintf(ex_s, "%g %g\n", x, real(E[j][0]));
		fprintf(ey_c, "%g %g\n", x, imag(E[j][1]));
		fprintf(ey_s, "%g %g\n", x, real(E[j][1]));
		fprintf(hx_c, "%g %g\n", x, imag(H[j][0]));
		fprintf(hx_s, "%g %g\n", x, real(H[j][0]));
		fprintf(hy_c, "%g %g\n", x, imag(H[j][1]));
		fprintf(hy_s, "%g %g\n", x, real(H[j][1]));	
		fprintf(hz_c, "%g %g\n", x, imag(H[j][2]));
		fprintf(hz_s, "%g %g\n", x, real(H[j][2]));
		fprintf(rok, "%g %g\n", x, rho[j]);
		fprintf(f_impedance, "%g %g\n", x, impedance[j]);
	}

	fclose(ex_c); fclose(ex_s); fclose(ey_c); fclose(ey_s);
	fclose(hx_c); fclose(hx_s); fclose(hy_c); fclose(hy_s); fclose(hz_c); fclose(hz_s);
	fclose(rok);
	fclose(f_impedance);
}

void Give_out_vec_mt::Write_B_to_files_for_harm_loop(int StartType,vector<int> &RecvPlsIgB,int ParamInd)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*B2d)[3];
	ifstream inf;
	ofstream outf;
	int ipls;
	char str[256];

	for(ipls=0;ipls<Subdomain::npls;ipls++)
	{
		if(!ParamInd)
			sprintf(str,"b3d_anom.%d",ipls+1);
		else
			sprintf(str,"d_b3d_anom.%d.%d",ipls+1,ParamInd);
		outf.open(str);
		outf<<setiosflags(ios_base::scientific)<<setprecision(14);
		for (i=0; i<d->n_pointresB; i++)
		{
			if(i<RecvPlsIgB[ipls] || i>=RecvPlsIgB[ipls+1]){continue;}

			outf<<MU_0*real(H[i+ipls*d->n_pointresB][0])<<' '
			<<MU_0*real(H[i+ipls*d->n_pointresB][1])<<' '
			<<MU_0*real(H[i+ipls*d->n_pointresB][2])<<' '
			<<MU_0*imag(H[i+ipls*d->n_pointresB][0])<<' '
			<<MU_0*imag(H[i+ipls*d->n_pointresB][1])<<' '
			<<MU_0*imag(H[i+ipls*d->n_pointresB][2])<<endl;
		}
		outf.close();
		outf.clear();
	}
}

void Give_out_vec_mt::Write_E_to_files_for_harm_loop(int StartType,vector<int> &RecvPlsIgE,int ParamInd)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*E2d)[3];
	ifstream inf;
	ofstream outf;
	int ipls;
	char str[256];

	for(ipls=0;ipls<Subdomain::npls;ipls++)
	{
		if(!ParamInd)
			sprintf(str,"e3d_anom.%d",ipls+1);
		else
			sprintf(str,"d_e3d_anom.%d.%d",ipls+1,ParamInd);
		outf.open(str);
		outf<<setiosflags(ios_base::scientific)<<setprecision(14);
		for (i=0; i<d->n_pointresE; i++)
		{
			if(i<RecvPlsIgE[ipls] || i>=RecvPlsIgE[ipls+1]){continue;}

			outf<<real(E[i+ipls*d->n_pointresE][0])<<' '
			<<real(E[i+ipls*d->n_pointresE][1])<<' '
			<<real(E[i+ipls*d->n_pointresE][2])<<' '
			<<imag(E[i+ipls*d->n_pointresE][0])<<' '
			<<imag(E[i+ipls*d->n_pointresE][1])<<' '
			<<imag(E[i+ipls*d->n_pointresE][2])<<endl;
		}
		outf.close();
		outf.clear();
	}
}

void Give_out_vec_mt::Write_result_to_files_for_harm_loop()
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	std::complex <double> (*B2d)[3];
	std::complex <double> (*E2d)[3];
	ifstream inf;
	ofstream outf;

	if((E2d = new std::complex<double>[d->n_pointresE][3])==0)
		Memory_allocation_error("E", "Write_result_to_files_for_harm_loop");

	inf.open("e2d");
	for (i=0; i<d->n_pointresE; i++)
	{
		inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
		E2d[i][0]=std::complex<double>(sin_x, cos_x); 
		E2d[i][1]=std::complex<double>(sin_y, cos_y); 
		E2d[i][2]=std::complex<double>(sin_z, cos_z); 
	}
	inf.close();
	inf.clear();

	if((B2d = new std::complex<double>[d->n_pointresB][3])==0)
		Memory_allocation_error("B", "Write_result_to_files_for_harm_loop");
		
	inf.open("b2d");
	for (i=0; i<d->n_pointresB; i++)
	{
		inf>>sin_x>>sin_y>>sin_z>>cos_x>>cos_y>>cos_z;
		B2d[i][0]=std::complex<double>(sin_x, cos_x); 
		B2d[i][1]=std::complex<double>(sin_y, cos_y); 
		B2d[i][2]=std::complex<double>(sin_z, cos_z); 
	}
	inf.close();
	inf.clear();	

	outf.open("e3d");
	for (i=0; i<d->n_pointresE; i++)
		outf<<real(E2d[i][0])+real(E[i][0])<<' '
			<<real(E2d[i][1])+real(E[i][1])<<' '
			<<real(E2d[i][2])+real(E[i][2])<<' '
			<<imag(E2d[i][0])+imag(E[i][0])<<' '
			<<imag(E2d[i][1])+imag(E[i][1])<<' '
			<<imag(E2d[i][2])+imag(E[i][2])<<endl; 
	outf.close();
	
	outf.open("b3d");
	for (i=0; i<d->n_pointresB; i++)
		outf<<real(B2d[i][0])+MU_0*real(H[i][0])<<' '
			<<real(B2d[i][1])+MU_0*real(H[i][1])<<' '
			<<real(B2d[i][2])+MU_0*real(H[i][2])<<' '
			<<imag(B2d[i][0])+MU_0*imag(H[i][0])<<' '
			<<imag(B2d[i][1])+MU_0*imag(H[i][1])<<' '
			<<imag(B2d[i][2])+MU_0*imag(H[i][2])<<endl; 
	outf.close();


	if(B2d) {delete [] B2d; B2d=NULL;}
	if(E2d) {delete [] E2d; E2d=NULL;}
}

int Give_out_vec_mt::Read_1d_field(char *fname_in, double &Ex_s, double &Ex_c, double &Ey_s, double &Ey_c, 
				  double &Hx_s, double &Hx_c, double &Hy_s, double &Hy_c)
{
	ifstream f_in;
	char buffer[10];

	f_in.open(fname_in);
	f_in >> buffer >> Ex_s;
	f_in >> buffer >> Ex_c;
	f_in >> buffer >> Ey_s;
	f_in >> buffer >> Ey_c;
	f_in >> buffer >> Hx_s;
	f_in >> buffer >> Hx_c;
	f_in >> buffer >> Hy_s;
	f_in >> buffer >> Hy_c;
	f_in.close();

	return 0;
}

int Give_out_vec_mt::GetNumberOfNodes()
{
	return d->kuzlov;
}
int Give_out_vec_mt::GetNumberOfElements()
{
	return d->kpar;
}
int Give_out_vec_mt::GetElementNodesNumber()
{
	return 8;
}
const pv::Point3D Give_out_vec_mt::GetNode(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyz[i_node][0];
	Point.y()=d->xyz[i_node][1];
	Point.z()=d->xyz[i_node][2];
	return Point;
}
const pv::Point3D Give_out_vec_mt::GetNodeTrue(const int& i_node)
{
	pv::Point3D Point;
	Point.x()=d->xyzt[i_node][0];
	Point.y()=d->xyzt[i_node][1];
	Point.z()=d->xyzt[i_node][2];
	return Point;
}
int Give_out_vec_mt::GetNodeNumberOnElement(const int& i_element, const int& i_node)
{
	return d->nver[i_element][i_node];
}
int Give_out_vec_mt::GetElementMaterial(const int& i_element)
{
	return d->nvkat[i_element];
}
int Give_out_vec_mt::GetTypeOfElement(const int& i_element)
{
	return d->nver[i_element][13];
}
double Give_out_vec_mt::GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)
{
	int i;
	double f[12];
	double x[8],y[8],z[8];
	double in[3],out[3],FieldOnElem;
	int isReal,isNoRot;

	FieldOnElem=0;

	isReal=r_type==vtAxSin ||r_type==vtAySin || r_type==vtAzSin || 
		r_type==vtRotxASin ||r_type==vtRotyASin || r_type==vtRotzASin;

	isNoRot=r_type==vtAxSin || r_type==vtAySin || r_type==vtAzSin || 
		r_type==vtAxCos ||r_type==vtAyCos || r_type==vtAzCos;
	
	for (i=0; i<12; i++){f[i] = v3dat[2*(tmap->ed[i_element][i])+!isReal+ipls_cur*tmap->n*2];}

	for (i=0; i<8; i++)
	{
		x[i] = d->xyz[d->nver[i_element][i]][0];
		y[i] = d->xyz[d->nver[i_element][i]][1];
		z[i] = d->xyz[d->nver[i_element][i]][2];
	}
	
	T_Brick L(x, y, z, d->nver[i_element][13]);

	in[0]=in[1]=in[2]=0.0;

	if(isNoRot)
	{
		L.Calc_value_inside_hex(f,in,out);

		if(r_type==vtAxSin || r_type==vtAxCos)FieldOnElem=out[0];
		if(r_type==vtAySin || r_type==vtAyCos)FieldOnElem=out[1];
		if(r_type==vtAzSin || r_type==vtAzCos)FieldOnElem=out[2];
	}
	else
	{
		L.Calc_rotor_inside_hex(f,in,out);

		if(r_type==vtRotxASin || r_type==vtRotxACos)FieldOnElem=out[0];
		if(r_type==vtRotyASin || r_type==vtRotyACos)FieldOnElem=out[1];
		if(r_type==vtRotzASin || r_type==vtRotzACos)FieldOnElem=out[2];
	}

	return FieldOnElem;
}
int Give_out_vec_mt::GetNumberOfResPoints(const Res3DValueType& r_type)
{
	return (r_type==vtWithDiscontinuity)? d->n_pointresE : d->n_pointresB;
}
pv::Point3D Give_out_vec_mt::GetResPoint(const Res3DValueType& r_type, const int& i_point)
{
	pv::Point3D Point;
	Point.x()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][0] : d->pointresB[i_point][0];
	Point.y()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][1] : d->pointresB[i_point][1];
	Point.z()=(r_type==vtWithDiscontinuity)? d->pointresE[i_point][2] : d->pointresB[i_point][2];
	return Point;
}

int * Give_out_vec_mt::GetPointerToRegular()
{
	return &(d->regular[0]);
}

int Give_out_vec_mt::GetXSize()
{
	return d->N_X;
}

int Give_out_vec_mt::GetYSize()
{
	return d->N_Y;
}

int Give_out_vec_mt::GetZSize()
{
	return d->N_Z;
}

double * Give_out_vec_mt::GetPointerToX()
{
	return &(d->Xcrd[0]);
}

double * Give_out_vec_mt::GetPointerToY()
{
	return &(d->Ycrd[0]);
}

double * Give_out_vec_mt::GetPointerToZ()
{
	return &(d->Zcrd[0]);
}

void Give_out_vec_mt::SaveResult(const Res3DValueType& r_type, const double& r_value, const int& j, const int& i_time, int ipls)
{
	double for_res_mu=MU_0;
	double for_res_w=d->nu*2.0*PI;

	if(r_type==vtRotxASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[ipls*d->n_pointresB+j][0]));
		H[ipls*d->n_pointresB+j][0]=ComplexValue;
	}
	else if(r_type==vtRotxACos){
		std::complex <double> ComplexValue(real(H[ipls*d->n_pointresB+j][0]),r_value/for_res_mu);
		H[ipls*d->n_pointresB+j][0]=ComplexValue;
	}
	else if(r_type==vtRotyASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[ipls*d->n_pointresB+j][1]));
		H[ipls*d->n_pointresB+j][1]=ComplexValue;
	}
	else if(r_type==vtRotyACos){
		std::complex <double> ComplexValue(real(H[ipls*d->n_pointresB+j][1]),r_value/for_res_mu);
		H[ipls*d->n_pointresB+j][1]=ComplexValue;
	}
	else if(r_type==vtRotzASin){
		std::complex <double> ComplexValue(r_value/for_res_mu,imag(H[ipls*d->n_pointresB+j][2]));
		H[ipls*d->n_pointresB+j][2]=ComplexValue;
	}
	else if(r_type==vtRotzACos){
		std::complex <double> ComplexValue(real(H[ipls*d->n_pointresB+j][2]),r_value/for_res_mu);
		H[ipls*d->n_pointresB+j][2]=ComplexValue;
	}
	else if(r_type==vtAxSin){
		std::complex <double> ComplexValue(real(E[ipls*d->n_pointresE+j][0]),-r_value*for_res_w);
		E[ipls*d->n_pointresE+j][0]=ComplexValue;
	}
	else if(r_type==vtAxCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[ipls*d->n_pointresE+j][0]));
		E[ipls*d->n_pointresE+j][0]=ComplexValue;
	}
	else if(r_type==vtAySin){
		std::complex <double> ComplexValue(real(E[ipls*d->n_pointresE+j][1]),-r_value*for_res_w);
		E[ipls*d->n_pointresE+j][1]=ComplexValue;
	}
	else if(r_type==vtAyCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[ipls*d->n_pointresE+j][1]));
		E[ipls*d->n_pointresE+j][1]=ComplexValue;
	}
	else if(r_type==vtAzSin){
		std::complex <double> ComplexValue(real(E[ipls*d->n_pointresE+j][2]),-r_value*for_res_w);
		E[ipls*d->n_pointresE+j][2]=ComplexValue;
	}
	else if(r_type==vtAzCos){
		std::complex <double> ComplexValue(r_value*for_res_w,imag(E[ipls*d->n_pointresE+j][2]));
		E[ipls*d->n_pointresE+j][2]=ComplexValue;
	}
	else
	{ 
		logfile<<"Unknown Res3DValueType in Give_out_vec_mt::SaveResult"<<endl;
	}
}

void Give_out_vec_mt::Write_NULL_B_to_files_for_harm_loop(int ParamInd)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	ifstream inf;
	ofstream outf;
	int ipls;
	char str[256];

	for(ipls=0;ipls<Subdomain::npls;ipls++)
	{
		if(!ParamInd)
			sprintf(str,"b3d_anom.%d",ipls+1);
		else
			sprintf(str,"d_b3d_anom.%d.%d",ipls+1,ParamInd);
		outf.open(str);
		outf.close();
		outf.clear();
	}
}

void Give_out_vec_mt::Write_NULL_E_to_files_for_harm_loop(int ParamInd)
{
	int i;
	double sin_x, sin_y, sin_z, cos_x, cos_y, cos_z; 
	ifstream inf;
	ofstream outf;
	int ipls;
	char str[256];

	for(ipls=0;ipls<Subdomain::npls;ipls++)
	{
		if(!ParamInd)
			sprintf(str,"e3d_anom.%d",ipls+1);
		else
			sprintf(str,"d_e3d_anom.%d.%d",ipls+1,ParamInd);
		outf.open(str);
		outf.close();
		outf.clear();
	}
}
