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

extern ofstream logfile;

MeshRZ::MeshRZ()
{
	l1 = NULL;
	sigma = NULL;
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
	if (mu) {delete [] mu; mu=NULL;}
	if (nvkat) {delete [] nvkat; nvkat=NULL;}
	if (type_of_rect) { delete [] type_of_rect; type_of_rect=NULL;}
	if (time) {delete [] time; time=NULL;}
	if (edges_2kr) {delete [] edges_2kr; edges_2kr=NULL;}
	if (reg) {delete [] reg; reg=NULL;}
	if (rm) {delete [] rm; rm=NULL;}
	if (zm) {delete [] zm; zm=NULL;}
}

int MeshRZ::ReadMeshRZ(char *path)
{
	char fname[128];
	In_Out r;
	int i, j;
	FILE *fp=NULL;
	double tmp1;
	int (*nver)[6]=NULL;

	sprintf(fname,"%s\\inf2tr.dat",path);
	r.Read_inf2tr(fname, &kuzlov, &krect, &kt1);

	n_materials = 0;
	sprintf(fname,"%s\\sigma",path);
	if((fp=fopen(fname, "r"))==0) Cannot_open_file(fname, "MeshRZ::ReadMeshRZ()");
	while(!feof(fp))
	{
		fscanf(fp, "%ld %lf ", &i, &tmp1);
		if(i > n_materials)
			n_materials = i;
	}
	fclose(fp);

	if((sigma = new double[n_materials])==0) Memory_allocation_error("sigma", "MeshRZ::ReadMeshRZ()");

	sprintf(fname,"%s\\sigma",path);
	if((fp=fopen(fname, "r"))==0) Cannot_open_file(fname, "MeshRZ::ReadMeshRZ()");
	for (i=0; i<n_materials; i++)
	{
		int n_of_current_material; 

		fscanf(fp, "%ld", &n_of_current_material);
		n_of_current_material--;
		fscanf(fp, "%lf", &sigma[n_of_current_material]);
	}
	fclose(fp);

	for(i=0;i<n_materials;i++)
	{
		if(sigma[i]<1e-6)
		{
			sigma[i]=1e-6;
		}
	}

	sigma_xy.resize(n_materials);
	for (i=0; i<n_materials; i++)
	{
		sigma_xy[i] = sigma[i];
	}

	sigma_z.resize(n_materials);

	sprintf(fname,"%s\\sigmaZ",path);
	if((fp=fopen(fname, "r"))==0) 
	{
		for (i=0; i<n_materials; i++)
		{
			sigma_z[i] = sigma[i];
		}
	}
	else
	{
		for (i=0; i<n_materials; i++)
		{
			int n_of_current_material; 

			fscanf(fp, "%ld", &n_of_current_material);
			n_of_current_material--;
			fscanf(fp, "%lf", &sigma_z[n_of_current_material]);
		}
		fclose(fp);

		for(i=0;i<n_materials;i++)
		{
			if(sigma_z[i]<1e-6)
			{
				sigma_z[i]=1e-6;
			}
		}
	}

	if((mu = new double[n_materials])==0) Memory_allocation_error("mu", "MeshRZ::ReadMeshRZ()");

	sprintf(fname,"%s\\mu",path);
	if((fp=fopen(fname, "r"))==0) Cannot_open_file(fname, "MeshRZ::ReadMeshRZ()");
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

	tsize2d = 0;
	n_nodes_c = kuzlov - tsize2d;

	if((l1 = new int[kt1])==0) Memory_allocation_error("l1", "MeshRZ::ReadMeshRZ()");

	sprintf(fname,"%s\\l1.dat",path);
	r.Read_Bin_File_Of_Long(fname, l1, kt1, 1);
	for(i=0; i<kt1; i++)
		l1[i]--;

	if((xy = new double[kuzlov][2])==0) Memory_allocation_error("rz", "MeshRZ::ReadMeshRZ()");
	sprintf(fname,"%s\\rz.dat",path);
	r.Read_Bin_File_Of_Double(fname, (double*)xy, kuzlov, 2);

	if((nver = new int[krect][6])==0) Memory_allocation_error("nver", "MeshRZ::ReadMeshRZ()");

	sprintf(fname,"%s\\nvtr.dat",path);
	r.Read_Bin_File_Of_Long(fname, (int*)nver, krect, 6);
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

	if((nvkat = new int[krect])==0) Memory_allocation_error("nvkat", "MeshRZ::ReadMeshRZ()");

	sprintf(fname,"%s\\nvkat2d.dat",path);
	r.Read_Bin_File_Of_Long(fname, nvkat, krect, 1);
	for(i=0; i<krect; i++)
		nvkat[i]--;

	ifstream inf;

	sprintf(fname,"%s\\rz.txt",path);
	inf.open(fname);
	if(!inf) return 1;
	inf>>qr>>qz;
	inf.close();
	inf.clear();

	if(!(rm=new double[qr])) return 1;
	if(!(zm=new double[qz])) return 1;

	sprintf(fname,"%s\\r.dat",path);
	inf.open(fname, ios::binary);
	if(!inf) return 1;
	for(i=0;i<qr;i++){inf>rm[i];}
	inf.close();
	inf.clear();

	sprintf(fname,"%s\\z.dat",path);
	inf.open(fname, ios::binary);
	if(!inf) return 1;
	for(i=0;i<qz;i++){inf>zm[i];}
	inf.close();
	inf.clear();

	nreg=krect;
	if(!(reg=new int[nreg])) return 1;
	for(i=0;i<nreg;i++){reg[i]=i+1;}

	return 0;
}

void MeshRZ::FixBoundCond()
{
	int i, k;
	double r;
	vector<bool> is_node_bound;

	is_node_bound.resize(kuzlov, false);

	for (i=0; i<kt1; i++)
		is_node_bound[l1[i]] = true;

	for (i=0; i<kuzlov; i++)
	{
		r = xy[i][0];
		if (r<=0.001)
		{
			is_node_bound[i] = false;
		}
	}

	cout << "kt1_old=" << kt1 << endl;

	kt1 = 0;
	for (i=0; i<kuzlov; i++)
	{
		if (is_node_bound[i])
			kt1++;
	}

	cout << "kt1_new=" << kt1 << endl;

	if (l1) {delete [] l1; l1=NULL;}
	l1 = new int[kt1];

	k = 0;
	for (i=0; i<kuzlov; i++)
	{
		if (is_node_bound[i])
		{
			l1[k] = i;
			k++;
		}
	}
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

	cout << "n_2kr=" << n_2kr << endl;

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
		
		if (fabs(r[0]-r0)<eps &&(z[0]>az || fabs(z[0]-az)<eps) &&(z[1]<bz || fabs(z[1]-bz)<eps))
		{
			edge[0] = v[0];
			edge[1] = v[1];
			edges.push_back(edge);									
		}		
	}

	n_2kr = (int)edges.size();
	
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


	bound_nodes.clear();
	bound_nodes.reserve(1000);

	bound_nodes.push_back(0);
	bound_nodes.push_back(1);
	bound_nodes.push_back(2);
	bound_nodes.push_back(3);
	bound_nodes.push_back(5);
	bound_nodes.push_back(6);
	bound_nodes.push_back(7);
	bound_nodes.push_back(8);

	kt1 = (int)bound_nodes.size();

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

int MeshRZ::FindNodeWithDeltaFunction(double z0, int *elem1, int *elem2)
{
	int i, j;
	int node, node_min;
	double r, r_min;
	double z, razn_z, razn_z_min;
	vector<int> nodes;
	const double eps=1e-5;
	int sz;
	int flag;
		 
	r_min = xy[0][0];

	for (i=1; i<kuzlov; i++)
	{
		r = xy[i][0];
		if (r < r_min)
		{
			r_min = r;
		}
	}

	for (i=0; i<kuzlov; i++)
	{
		r = xy[i][0];
		if (fabs(r_min - r)<eps)
		{
			nodes.push_back(i);
		}
	}

	sz = (int)nodes.size();


	node_min = nodes[0];
	z = xy[node_min][1];
	razn_z_min = fabs(z - z0);

	for (i=1; i<sz; i++)
	{
		node = nodes[i];
		z = xy[node][1];
		razn_z = fabs(z - z0);
		if (razn_z < razn_z_min)
		{
			node_min = node;
			razn_z_min = razn_z;
		}
	}

	flag = 0;

	for (i=0; i<krect; i++)
	{
		if (flag == 2)
			break;

		for (j=0; j<4; j++)
		{
			node = nvtr[i][j];

			if (node == node_min)
			{
				if (flag == 0)
				{
					*elem1 = i;
					flag++;
					break;
				} 
				else if (flag == 1)
				{
					*elem2 = i;
					flag++;
					break;
				}
			}
		}
	}

	return node_min;
}
