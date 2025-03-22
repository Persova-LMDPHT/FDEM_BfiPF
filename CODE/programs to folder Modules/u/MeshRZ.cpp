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
 *  Written by  D.Sc. Denis V. Vagin                                           
 *  Novosibirsk State Technical University,                                    
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                          
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)           
 *  Version 2.0 December 10, 2024                                              
*/

#include "stdafx.h"
#include "MeshRZ.h"
#include "in_out.h"
#include "For_Solvers.h"
extern ofstream logfile;

MeshRZ::MeshRZ()
{
	l1 = NULL;
	sigma = NULL;
	sigmaZ = NULL;
	dpr = NULL;
	mu = NULL;
	xy = NULL;
	nvtr = NULL;
	type_of_rect = NULL;
	nvkat = NULL;
	time = NULL;
	edges_2kr = NULL;
	reg = NULL;
	rm = NULL;
	zm = NULL;
}

MeshRZ::~MeshRZ()
{
	if (nvtr) {delete [] nvtr; nvtr=NULL;}
	if (xy) {delete [] xy; xy=NULL;}
	if (l1) {delete [] l1; l1=NULL;}
	if (sigma) {delete [] sigma; sigma=NULL;}
	if (sigmaZ) {delete [] sigmaZ; sigmaZ=NULL;}
	if (dpr) {delete [] dpr; dpr=NULL;}
	if (mu) {delete [] mu; mu=NULL;}
	if (nvkat) {delete [] nvkat; nvkat=NULL;}

	if (type_of_rect) { delete [] type_of_rect; type_of_rect=NULL;}

	if (time) {delete [] time; time=NULL;}

	if (edges_2kr) {delete [] edges_2kr; edges_2kr=NULL;}

	if (reg) {delete [] reg; reg=NULL;}
	if (rm) {delete [] rm; rm=NULL;}
	if (zm) {delete [] zm; zm=NULL;}
}

int MeshRZ::ReadMeshRZ()
{
	In_Out r;
	int i, j;
	FILE *fp=NULL;
	double tmp1;
	int (*nver)[6]=NULL;

	r.Read_inf2tr("Ax\\inf2tr.dat", &kuzlov, &krect, &kt1);
	n_materials = 0;
	if((fp=fopen("Ax\\sigma", "r"))==0) Cannot_open_file("Ax\\sigma", "MeshRZ::ReadMeshRZ()");
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf ", &i, &tmp1);
		if(i > n_materials)
			n_materials = i;
	}
	fclose(fp);

	if((sigma = new double[n_materials])==0) Memory_allocation_error("sigma", "MeshRZ::ReadMeshRZ()");

	if((fp=fopen("Ax\\sigma", "r"))==0) Cannot_open_file("Ax\\sigma", "MeshRZ::ReadMeshRZ()");
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &sigma[n_of_current_material]);
		if(sigma[n_of_current_material]<1e-8){sigma[n_of_current_material]=1e-8;}
	}
	fclose(fp);

	if((sigmaZ = new double[n_materials])==0) Memory_allocation_error("sigmaZ", "MeshRZ::ReadMeshRZ()");

	if((fp=fopen("Ax\\sigmaZ", "r"))==0) Cannot_open_file("Ax\\sigmaZ", "MeshRZ::ReadMeshRZ()");
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &sigmaZ[n_of_current_material]);
		if(sigmaZ[n_of_current_material]<1e-8){sigmaZ[n_of_current_material]=1e-8;}
	}
	fclose(fp);

	if((dpr = new double[n_materials])==0) Memory_allocation_error("dpr", "MeshRZ::ReadMeshRZ()");

	if((fp=fopen("Ax\\dpr", "r"))==0) Cannot_open_file("Ax\\dpr", "MeshRZ::ReadMeshRZ()");
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &dpr[n_of_current_material]);
		dpr[n_of_current_material]*=8.84194128288307421e-12;
	}
	fclose(fp);

	if((mu = new double[n_materials])==0) Memory_allocation_error("mu", "MeshRZ::ReadMeshRZ()");

	if((fp=fopen("Ax\\mu", "r"))==0) Cannot_open_file("Ax\\mu", "MeshRZ::ReadMeshRZ()");
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material; 
		int retcode;

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &mu[n_of_current_material]);
		mu[n_of_current_material] *= MU_0;
	}
	fclose(fp);

	r.Read_int_From_Txt_File("Ax\\tsize.dat", &tsize2d);

	n_nodes_c = kuzlov - tsize2d;

	if((l1 = new int[kt1])==0) Memory_allocation_error("l1", "MeshRZ::ReadMeshRZ()");

	r.Read_Bin_File_Of_Long("Ax\\l1.dat", l1, kt1, 1);
	for(i=0; i<kt1; i++)
		l1[i]--;

	if((xy = new double[kuzlov][2])==0) Memory_allocation_error("rz", "MeshRZ::ReadMeshRZ()");
	r.Read_Bin_File_Of_Double("Ax\\rz.dat", (double*)xy, kuzlov, 2);

	if((nver = new int[krect][6])==0) Memory_allocation_error("nver", "MeshRZ::ReadMeshRZ()");

	r.Read_Bin_File_Of_Long("Ax\\nvtr.dat", (int*)nver, krect, 6);
	for(i=0; i<krect; i++)
		for(j=0; j<6; j++)
			nver[i][j]--;

	if((type_of_rect = new int[krect])==0) Memory_allocation_error("type_of_rect", "MeshRZ::ReadMeshRZ()");

	for(i=0; i<krect; i++)
	{
		type_of_rect[i] = Rectangle_or_quadrilateral(i, xy, nver);
	}

	if((nvtr = new int[krect][4])==0) Memory_allocation_error("nvtr", "MeshRZ::ReadMeshRZ()");

	for(i=0; i<krect; i++)
	{
		nvtr[i][0] = nver[i][2];
		nvtr[i][1] = nver[i][3];
		nvtr[i][2] = nver[i][0];
		nvtr[i][3] = nver[i][1];
	}

	if(nver) {delete [] nver; nver=NULL;}

	if((nvkat = new int[krect])==0) Memory_allocation_error("nvkat2d", "MeshRZ::ReadMeshRZ()");

	r.Read_Bin_File_Of_Long("Ax\\nvkat2d.dat", nvkat, krect, 1);
	for(i=0; i<krect; i++)
		nvkat[i]--;

	ifstream inf;

	inf.open("Ax\\rz.txt");
	if (!inf) return 1;
	inf>>qr>>qz;
	inf.close();
	inf.clear();

	nreg=(qr-1)*(qz-1);

	if(!(reg=new int[nreg])) return 1;
	if(!(rm=new double[qr])) return 1;
	if(!(zm=new double[qz])) return 1;

	inf.open("Ax\\r.dat", ios::binary);
	if (!inf) return 1;
	for(i=0;i<qr;i++){inf>rm[i];}
	inf.close();
	inf.clear();

	inf.open("Ax\\z.dat", ios::binary);
	if (!inf) return 1;
	for(i=0;i<qz;i++){inf>zm[i];}
	inf.close();
	inf.clear();

	if(tsize2d>0)
	{
		inf.open("Ax\\reg.dat", ios::binary);
		if (!inf) return 1;
		for(i=0;i<nreg;i++){inf>reg[i];}
		inf.close();
		inf.clear();
	}
	else
	{
		for(i=0;i<nreg;i++){reg[i]=i+1;}
	}

	ofstream fout("rm");
	fout << qr << endl;
	for(i=0;i<qr;i++)
	{
		fout << rm[i] << endl;
	}
	fout.close();
	fout.clear();

	fout.open("zm");
	fout << qz << endl;
	for(i=0; i<qz; i++)
	{
		fout << zm[i] << endl;
	}
	fout.close();
	fout.clear();

	return 0;
}

int MeshRZ::Rectangle_or_quadrilateral(int num, double (*xy)[2], int (*nver)[6])
{
	double x[4], y[4];
	const double eps = 1e-8;
	int i;
	int type;

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

int MeshRZ::Read2Kr(char *f_n, char *f_2kr)
{
	n_2kr = 0;
	In_Out r;
	int i;

	r.Read_int_From_Txt_File(f_n, &n_2kr);

	if (edges_2kr) {delete [] edges_2kr; edges_2kr=NULL;}
	edges_2kr = new int[n_2kr][2];

	r.Read_Txt_File_Of_Long(f_2kr, (int*)edges_2kr, n_2kr, 2);

	for (i=0; i<n_2kr; i++)
	{
		edges_2kr[i][0]--;
		edges_2kr[i][1]--;
	}

	return 0;
}

int MeshRZ::Find2kr(double az, double bz, double r0, double eps)
{
	int i, j;
	int v[2];
	double r[2];
	double z[2];
	vector< vector<int> > edges;
	vector<int> edge;

	vector<int> nodes_tmp;

	edges.clear();
	edges.reserve(100);

	edge.resize(2);

	if (az > bz)
	{
		double temp;
		temp = az;
		az = bz;
		bz =temp;
	}
	
	for (i=0; i<krect; i++)
	{
		v[0] = nvtr[i][0];	
		v[1] = nvtr[i][2];

		for (j=0; j<2; j++)
		{
			r[j] = xy[v[j]][0];
			z[j] = xy[v[j]][1];
		}
		
		if (fabs(r[0]-r0)<eps && (z[0]>az || fabs(z[0]-az)<eps) && (z[1]<bz || fabs(z[1]-bz)<eps))
		{
			edge[0] = v[0];
			edge[1] = v[1];
			edges.push_back(edge);									
		}		
	}

	n_2kr = edges.size();
	
	if (edges_2kr) {delete [] edges_2kr; edges_2kr = NULL;}
	edges_2kr = new int[n_2kr][2];

	for (i=0; i<n_2kr; i++)
	{
		edges_2kr[i][0] = edges[i][0];
		edges_2kr[i][1] = edges[i][1];
	}

	nodes_tmp.resize(n_2kr*2);
	for (i=0; i<n_2kr; i++)
	{
		for (j=0; j<2; j++)
		{
			nodes_tmp[i*2+j] = edges[i][j];
		}
	}

	nodes_2kr.clear();
	std::sort(nodes_tmp.begin(), nodes_tmp.end());
	std::unique_copy(nodes_tmp.begin(), nodes_tmp.end(), back_inserter(nodes_2kr));

	In_Out R;

	R.Write_int_To_Txt_File("n_2kr.txt", n_2kr);
	R.Write_Txt_File_Of_Long("edges_2kr.txt", (int*)edges_2kr, n_2kr, 2);

	return 0;
}

int MeshRZ::Find1kr()
{
	int i;
	double r_min=rm[0], r_max=rm[qr-1], z_min=zm[0], z_max=zm[qz-1];
	double r;
	double z;
	vector <int> bound_nodes;
	double eps=1e-6;

	for (i=0; i<kuzlov; i++)
	{
		r = xy[i][0];
		z = xy[i][1];

		if (fabs(r-r_min)<eps || fabs(r-r_max)<eps || fabs(z-z_min)<eps || fabs(z-z_max)<eps)
		{
			{
				bound_nodes.push_back(i);
			}
		}
	}

	kt1 = bound_nodes.size();

	if (l1) {delete [] l1; l1=NULL;}
	l1 = new int[kt1];

	for (i=0; i<kt1; i++)
	{
		l1[i] = bound_nodes[i];
	}

	In_Out R;
	R.Write_Txt_File_Of_Long("l1.txt", l1, kt1, 1);

	{
		ofstream fout;
		fout.open("l1_coords.txt");

		for (i=0; i<kt1; i++)
		{
			fout << l1[i] << '\t' << xy[l1[i]][0] << '\t' << xy[l1[i]][1] << '\n';
		}
	}
	return 0;
}

int MeshRZ::Find1kr(double eps)
{
	int i;
	double r_min, r_max, z_min, z_max;
	double r;
	double z;
	vector <int> bound_nodes;

	r_min = r_max = xy[0][0];
	z_min = z_max = xy[0][1];

	for (i=0; i<kuzlov; i++)
	{
		r = xy[i][0];
		z = xy[i][1];

		if (r < r_min)
			r_min = r;

		if (r > r_max)
			r_max = r;

		if (z < z_min)
			z_min = z;

		if (z > z_max)
			z_max = z;
	}

	for (i=0; i<kuzlov; i++)
	{
		r = xy[i][0];
		z = xy[i][1];

		if (fabs(r-r_min)<eps || fabs(r-r_max)<eps || fabs(z-z_min)<eps || fabs(z-z_max)<eps)
		{
			{
				bound_nodes.push_back(i);
			}
		}
	}


	kt1 = bound_nodes.size();

	if (l1) {delete [] l1; l1=NULL;}
	l1 = new int[kt1];

	for (i=0; i<kt1; i++)
	{
		l1[i] = bound_nodes[i];
	}

	In_Out R;
	R.Write_Txt_File_Of_Long("l1.txt", l1, kt1, 1);

	{
		ofstream fout;
		fout.open("l1_coords.txt");

		for (i=0; i<kt1; i++)
		{
			fout << l1[i] << '\t' << xy[l1[i]][0] << '\t' << xy[l1[i]][1] << '\n';
		}
	}


	return 0;
}

int MeshRZ::Read_Az_Bz(char *fname)
{
	ifstream fin;
	double tmp;

	fin.open(fname);

	if (!fin)
		Cannot_open_file(fname, "MeshRZ::Read_Az_Bz");

	fin >> tmp;
	fin >> tmp;
	fin >> az_line;
	fin >> tmp;
	fin >> tmp;
	fin >> bz_line;

	fin.close();

	return 0;
}
