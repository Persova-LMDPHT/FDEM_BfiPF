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
 *  This file contains a code functions for finite element mesh storing
 *  
 *  Written by D.Sc. Denis V. Vagin     
 *  Novosibirsk State Technical University,                                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                        
 *  Version 2.0 December 10, 2024                                                             
*/

#include "stdafx.h"
#include "rect_postproc_2d.h"
#include "in_out.h"
#include "For_Solvers.h"
#include "rect_local_matrix_2d.h"
#include "iobinary.h"

extern ofstream logfile;

Rect_preproc_data::Rect_preproc_data()
{
	mu2d=NULL;
	mu0=NULL;
	sigma2d=NULL;
	sigma0=NULL;
	l1=NULL;
	nvk1=NULL;
	nvtr=NULL;
	nvkat=NULL;
	xy=NULL;
	type_of_rect=NULL;

	coords_1d=NULL;
	sin_1d=NULL;
	cos_1d=NULL;

	pointres=NULL;
	xyzVectorE0=NULL;
	xyzVectorB=NULL;
	xyzVectorB2d=NULL;
	xyzVectorE2d=NULL;

	mtr3d2d=NULL;

	n_pointres = 0;
	n_xyzVectorE0 = 0;
	n_xyzVectorB = 0;
	n_xyzVectorB2d = 0;
	n_xyzVectorE2d = 0;

	pType=ptNone;
}

Rect_preproc_data::~Rect_preproc_data()
{
	if(mu2d) { delete [] mu2d; mu2d=NULL;}
	if(mu0) { delete [] mu0; mu0=NULL;}
	if(sigma2d) { delete [] sigma2d; sigma2d=NULL;}
	if(sigma0) { delete [] sigma0; sigma0=NULL;}
	if(l1) { delete [] l1; l1=NULL;}
	if(nvk1) { delete [] nvk1; nvk1=NULL;}

	if(nvtr) { delete [] nvtr; nvtr=NULL;}
	if(nvkat) { delete [] nvkat; nvkat=NULL;}
	if(xy) { delete [] xy; xy=NULL;}
	if(type_of_rect) { delete [] type_of_rect; type_of_rect=NULL;}

	if(coords_1d) { delete [] coords_1d; coords_1d=NULL;}
	if(sin_1d) { delete [] sin_1d; sin_1d=NULL;}
	if(cos_1d) { delete [] cos_1d; cos_1d=NULL;}

	if(pointres) { delete [] pointres; pointres=NULL;}
	if(xyzVectorE0) { delete [] xyzVectorE0; xyzVectorE0=NULL;}
	if(xyzVectorB) { delete [] xyzVectorB; xyzVectorB=NULL;}
	if(xyzVectorB2d) { delete [] xyzVectorB2d; xyzVectorB2d=NULL;}
	if(xyzVectorE2d) { delete [] xyzVectorE2d; xyzVectorE2d=NULL;}

	if (mtr3d2d) { delete [] mtr3d2d; mtr3d2d=NULL; }
}

int Rect_preproc_data::Read_1d()
{
	In_Out R;

	R.Read_Long_From_Txt_File("param", &n_1d);

	if((coords_1d = new double[n_1d])==0) Memory_allocation_error("coords_1d", "Rect_preproc_data::Read_1d");
	if((sin_1d = new double[n_1d])==0) Memory_allocation_error("sin_1d", "Rect_preproc_data::Read_1d");;
	if((cos_1d = new double[n_1d])==0) Memory_allocation_error("cos_1d", "Rect_preproc_data::Read_1d");;

	R.Read_1d_data(n_1d, coords_1d, sin_1d, cos_1d);

	return 0;
}

int Rect_preproc_data::Read_pointres()
{
	FILE *fp=NULL;
	long i;
	double tmp1, tmp2;

	if((fp=fopen("pointres", "r"))==0) Cannot_open_file("pointres", "Rect_preproc_data::Read_pointres");
	
	fscanf(fp, "%ld", &n_pointres);

	if((pointres = new double[n_pointres])==0) Memory_allocation_error("pointres","Rect_preproc_data::Read_pointres");

	for(i=0; i<n_pointres; i++)
        fscanf(fp, "%lf %lf %lf", &pointres[i], &tmp1, &tmp2);

	fclose(fp);
	return 0;
}

int Rect_preproc_data::Read_rect_preproc_data()
{
	In_Out R;
	FILE *fp=NULL;
	long i, j;
	double alph;
	double tmp1, tmp2;
	long (*nver)[6]=NULL;

	R.Read_Double_From_Txt_File("alfa", &alph);
	if(fabs(alph)<1e-3)
	{
		alfa = 0;
	}
	else
	{
		alfa = 1;
	}

	if((fp=fopen("sig2d", "r"))==0) Cannot_open_file("sig2d", "Rect_preproc_data::Read_rect_preproc_data");

	fscanf(fp, "%ld", &i);
	n_materials = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf %lf", &i, &tmp1, &tmp2);
			if(i > n_materials)
				n_materials = i;
	}
	fclose(fp);

	if((fp=fopen("sig2d", "r"))==0) Cannot_open_file("sig2d", "Rect_preproc_data::Read_rect_preproc_data");

	fscanf(fp, "%ld", &i);

	if((sigma2d = new double[n_materials])==0) Memory_allocation_error("sigma2d", "Rect_preproc_data::Read_rect_preproc_data");
	if((sigma0 = new double[n_materials])==0) Memory_allocation_error("sigma0", "Rect_preproc_data::Read_rect_preproc_data");

	for(i=0; i<n_materials; i++)
	{
		long n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf %lf", &sigma2d[n_of_current_material], &sigma0[n_of_current_material]);
	}
	fclose(fp);

	if((mu2d = new double[n_materials])==0) Memory_allocation_error("mu2d", "Rect_preproc_data::Read_rect_preproc_data");
	if((mu0 = new double[n_materials])==0) Memory_allocation_error("mu0", "Rect_preproc_data::Read_rect_preproc_data");

	for(i=0; i<n_materials; i++)
	{
		mu2d[i] = 4.0*3.1415926535*1E-7;
		mu0[i] = 4.0*3.1415926535*1E-7;
	}
	
	R.Read_Double_From_Txt_File("nu", &nu);
	R.Read_inf2tr("inf2tr.dat", &kuzlov, &krect, &kt1);
	R.Read_Long_From_Txt_File("tsize.dat", &tsize2d);

	n_nodes_c = kuzlov - tsize2d;

	if((l1 = new long[kt1])==0) Memory_allocation_error("l1", "Rect_preproc_data::Read_rect_preproc_data");

	R.Read_Bin_File_Of_Long("l2d.dat", l1, kt1, 1);
	for(i=0; i<kt1; i++)
		l1[i]--;

	if((xy = new double[kuzlov][2])==0) Memory_allocation_error("xy", "Rect_preproc_data::Read_rect_preproc_data");

	R.Read_Bin_File_Of_Double("xz.dat", (double*)xy, kuzlov, 2);

	if((nver = new long[krect][6])==0) Memory_allocation_error("nver", "Rect_preproc_data::Read_rect_preproc_data");

	R.Read_Bin_File_Of_Long("nvtr.dat", (long*)nver, krect, 6);
	for(i=0; i<krect; i++)
		for(j=0; j<6; j++)
			nver[i][j]--;

	if((type_of_rect = new long[krect])==0) Memory_allocation_error("type_of_rect", "Rect_preproc_data::Read_rect_preproc_data");

	for(i=0; i<krect; i++)
	{
		type_of_rect[i] = Rectangle_or_quadrilateral(i, xy, nver);
	}

	if((nvtr = new long[krect][4])==0) Memory_allocation_error("nvtr", "Rect_preproc_data::Read_rect_preproc_data");

	for(i=0; i<krect; i++)
		for(j=0; j<4; j++)
			nvtr[i][j] = nver[i][j];

	if(nver) {delete [] nver; nver=NULL;}

	if((nvkat = new long[krect])==0) Memory_allocation_error("nvkat", "Rect_preproc_data::Read_rect_preproc_data");

	R.Read_Bin_File_Of_Long("nvkat.dat", nvkat, krect, 1);
	for(i=0; i<krect; i++)
		nvkat[i]--;

	return 0;
}

long Rect_preproc_data::Rectangle_or_quadrilateral(long num, double (*xy)[2], long (*nver)[6])
{
	double x[4], y[4];
	const double eps = 1e-2;
	long i;
	long type;

	type = nver[num][5];

	for(i=0; i<4; i++)
	{
		x[i] = xy[nver[num][i]][0];
		y[i] = xy[nver[num][i]][1];
	}

	if(fabs(y[0]-y[1])>eps) return type+5; 
	if(fabs(y[2]-y[3])>eps) return type+5; 
	if(fabs(x[0]-x[2])>eps) return type+5; 
	if(fabs(x[1]-x[3])>eps) return type+5; 

	return type;
}

int Rect_preproc_data::Read_mesh_rz()
{
	In_Out R;
	FILE *fp=NULL;
	ifstream fin;
	long i, j, k;
	long (*nver)[6]=NULL;

	long max_material, tmp;
	double temp1;

	if((fp=fopen("sigma", "r"))==0)
	{
		Cannot_open_file("sigma", "Hex_preproc_data::Read_mesh_rz");
		return 1;
	}
	max_material = 0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf", &tmp, &temp1);
		if(tmp > max_material)
			max_material = tmp;
	}
	n_materials = max_material;
	fclose(fp);

	if((fp=fopen("sigma", "r"))==0) Cannot_open_file("sigma", "Rect_preproc_data::Read_mesh_rz");

	if((sigma2d = new double[n_materials])==0) Memory_allocation_error("sigma2d", "Rect_preproc_data::Read_mesh_rz");

	for(i=0; i<n_materials; i++)
	{
		long n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &sigma2d[n_of_current_material]);
	}
	fclose(fp);

	if((fp=fopen("dpr", "r"))==0) Cannot_open_file("dpr", "Rect_preproc_data::Read_mesh_rz");

	if((dpr2d = new double[n_materials])==0) Memory_allocation_error("dpr2d", "Rect_preproc_data::Read_mesh_rz");
	
	for(i=0; i<n_materials; i++)
	{
		long n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &dpr2d[n_of_current_material]);
	}
	fclose(fp);
	
	if((mu2d = new double[n_materials])==0) Memory_allocation_error("mu2d", "Rect_preproc_data::Read_mesh_rz");
	if((mu0 = new double[n_materials])==0) Memory_allocation_error("mu0", "Rect_preproc_data::Read_mesh_rz");

	for(i=0; i<n_materials; i++)
	{
		mu2d[i] = 4.0*3.1415926535*1E-7;
		mu0[i] = 4.0*3.1415926535*1E-7;
	}

	if((fp=fopen("mu", "r"))==0) Cannot_open_file("mu", "Rect_preproc_data::Read_mesh_rz");
	
	for(i=0; i<n_materials; i++)
	{
		long n_of_current_material; 
		double _mu2d;

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &_mu2d);
		mu2d[n_of_current_material]*=_mu2d;
	}
	fclose(fp);

	R.Read_Double_From_Txt_File("nu", &nu);
	R.Read_inf2tr("inf2tr.dat", &kuzlov, &krect, &kt1);
	R.Read_Long_From_Txt_File("tsize.dat", &tsize2d);
	n_nodes_c = kuzlov - tsize2d;

	if((l1 = new long[kt1])==0) Memory_allocation_error("l1", "Rect_preproc_data::Read_mesh_rz");

	R.Read_Bin_File_Of_Long("l1.dat", l1, kt1, 1);
	for(i=0; i<kt1; i++)
		l1[i]--;

	if (pType==ptHphi)
	{
		if((nvk1 = new short[kt1])==0) Memory_allocation_error("nvk1", "Rect_preproc_data::Read_mesh_rz");
		R.Read_Bin_File_Of_Short("nvk1.dat", nvk1, kt1, 1);
	}

	if((xy = new double[kuzlov][2])==0) Memory_allocation_error("xy", "Rect_preproc_data::Read_mesh_rz");

	R.Read_Bin_File_Of_Double("rz.dat", (double*)xy, kuzlov, 2);

	if((nver = new long[krect][6])==0) Memory_allocation_error("nver", "Rect_preproc_data::Read_mesh_rz");

	R.Read_Bin_File_Of_Long("nvtr.dat", (long*)nver, krect, 6);
	for(i=0; i<krect; i++)
		for(j=0; j<6; j++)
			nver[i][j]--;

	if((type_of_rect = new long[krect])==0) Memory_allocation_error("type_of_rect", "Rect_preproc_data::Read_mesh_rz");

	for(i=0; i<krect; i++)
	{
		type_of_rect[i] = Rectangle_or_quadrilateral(i, xy, nver);
	}

	if((nvtr = new long[krect][4])==0) Memory_allocation_error("nvtr", "Rect_preproc_data::Read_mesh_rz");

	for(i=0; i<krect; i++)
		for(j=0; j<4; j++)
			nvtr[i][j] = nver[i][j];

	if(nver) {delete [] nver; nver=NULL;}

	if((nvkat = new long[krect])==0) Memory_allocation_error("nvkat", "Rect_preproc_data::Read_mesh_rz");

	R.Read_Bin_File_Of_Long("nvkat2d.dat", nvkat, krect, 1);
	for(i=0; i<krect; i++)
		nvkat[i]--;

	fin.open("rz.txt");
	if(!fin) Cannot_open_file("rz.txt", "Rect_preproc_data::Read_mesh_rz");
	fin>>qr>>qz;
	fin.close();
	fin.clear();

	rm.resize(qr);
	zm.resize(qz);

	fin.open("r.dat", ios::binary);
	if(!fin) Cannot_open_file("r.dat", "Rect_preproc_data::Read_mesh_rz");
	for(i=0;i<qr;i++){fin>rm[i];}
	fin.close();
	fin.clear();

	fin.open("z.dat", ios::binary);
	if(!fin) Cannot_open_file("z.dat", "Rect_preproc_data::Read_mesh_rz");
	for(i=0;i<qz;i++){fin>zm[i];}
	fin.close();
	fin.clear();

	nreg=krect;
	reg.resize(nreg);
	for(i=0;i<nreg;i++){reg[i]=i+1;}

	if((fp=fopen("mtr3d2d", "r"))==0)
	{
		Cannot_open_file("mtr3d2d", "Hex_preproc_data::Read_mesh_rz");
		return 1;
	}
	k=0;
	while(!feof(fp))
	{
		fscanf(fp, "%ld %ld", &tmp, &j);
		if(tmp>k)
			k = tmp;
	}
	fclose(fp);
	
	fin.open("mtr3d2d");
	if (fin)
	{
		mtr3d2d = new int[k];
		for (i=0; i<k; i++)
		{
			fin>>tmp>>j;
			mtr3d2d[tmp-1]=j;
		}
		fin.close();
	}
	
	return 0;
}

int Rect_preproc_data::Read_xyzVectorE0()
{
	FILE *fp=NULL;
	long i;

	if((fp=fopen("xyzVectorE0", "r"))==0)
	{
		Cannot_open_file_but_continue("xyzVectorE0", "Rect_preproc_data::Read_xyzVectorE0");
		return -1;
	}

	fscanf_s(fp, "%ld", &n_xyzVectorE0);

	if(xyzVectorE0) {delete [] xyzVectorE0; xyzVectorE0=NULL;}
	if((xyzVectorE0 = new double[n_xyzVectorE0][3])==0) Memory_allocation_error("xyzVectorE0","Rect_preproc_data::Read_xyzVectorE0");

	for(i=0; i<n_xyzVectorE0; i++)
		fscanf_s(fp, "%lf %lf %lf", &xyzVectorE0[i][0], &xyzVectorE0[i][1], &xyzVectorE0[i][2]);

	fclose(fp);
	return 0;
}

int Rect_preproc_data::Read_xyzVectorB()
{
	FILE *fp=NULL;
	long i;

	if((fp=fopen("xyzVectorB", "r"))==0)
	{
		Cannot_open_file_but_continue("xyzVectorB", "Rect_preproc_data::Read_xyzVectorB");
		return -1;
	}

	fscanf_s(fp, "%ld", &n_xyzVectorB);

	if(xyzVectorB) {delete [] xyzVectorB; xyzVectorB=NULL;}
	if((xyzVectorB = new double[n_xyzVectorB][3])==0) Memory_allocation_error("xyzVectorB","Rect_preproc_data::Read_xyzVectorB");

	for(i=0; i<n_xyzVectorB; i++)
		fscanf_s(fp, "%lf %lf %lf", &xyzVectorB[i][0], &xyzVectorB[i][1], &xyzVectorB[i][2]);

	fclose(fp);
	return 0;
}

int Rect_preproc_data::Read_xyzVectorE()
{
	FILE *fp=NULL;
	long i;

	if((fp=fopen("xyzVectorE", "r"))==0)
	{
		Cannot_open_file_but_continue("xyzVectorE", "Rect_preproc_data::Read_xyzVectorE");
		return -1;
	}

	fscanf_s(fp, "%ld", &n_xyzVectorE);

	if(xyzVectorE) {delete [] xyzVectorE; xyzVectorE=NULL;}
	if((xyzVectorE = new double[n_xyzVectorE][3])==0) Memory_allocation_error("xyzVectorE","Rect_preproc_data::Read_xyzVectorE");

	for(i=0; i<n_xyzVectorE; i++)
		fscanf_s(fp, "%lf %lf %lf", &xyzVectorE[i][0], &xyzVectorE[i][1], &xyzVectorE[i][2]);

	fclose(fp);
	return 0;
}

int Rect_preproc_data::Read_xyzVectorB2d()
{
	FILE *fp=NULL;
	long i;

	if((fp=fopen("xyzVectorB2d", "r"))==0)
	{
		Cannot_open_file_but_continue("xyzVectorB2d", "Rect_preproc_data::Read_xyzVectorB2d");
		return -1;
	}

	fscanf_s(fp, "%ld", &n_xyzVectorB2d);

	if(xyzVectorB2d) {delete [] xyzVectorB2d; xyzVectorB2d=NULL;}
	if((xyzVectorB2d = new double[n_xyzVectorB2d][3])==0) Memory_allocation_error("xyzVectorB2d","Rect_preproc_data::Read_xyzVectorB2d");

	for(i=0; i<n_xyzVectorB2d; i++)
		fscanf_s(fp, "%lf %lf %lf", &xyzVectorB2d[i][0], &xyzVectorB2d[i][1], &xyzVectorB2d[i][2]);

	fclose(fp);
	return 0;
}

int Rect_preproc_data::Read_xyzVectorE2d()
{
	FILE *fp=NULL;
	long i;

	if((fp=fopen("xyzVectorE2d", "r"))==0)
	{
		Cannot_open_file_but_continue("xyzVectorE2d", "Rect_preproc_data::Read_xyzVectorE2d");
		return -1;
	}

	fscanf_s(fp, "%ld", &n_xyzVectorE2d);

	if(xyzVectorE2d) {delete [] xyzVectorE2d; xyzVectorE2d=NULL;}
	if((xyzVectorE2d = new double[n_xyzVectorE2d][3])==0) Memory_allocation_error("xyzVectorE2d","Rect_preproc_data::Read_xyzVectorE2d");

	for(i=0; i<n_xyzVectorE2d; i++)
		fscanf_s(fp, "%lf %lf %lf", &xyzVectorE2d[i][0], &xyzVectorE2d[i][1], &xyzVectorE2d[i][2]);

	fclose(fp);
	return 0;
}
