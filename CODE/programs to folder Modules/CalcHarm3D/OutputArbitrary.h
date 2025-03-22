/*
 * GENERAL REMARKS
 * 
 *  This code is freely available under the following conditions:
 * 
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *  
 *                          
 *  This file contains headers for some basic routines for output of values of the desired function
 *  at an arbitrary point in the computational domain
 * 
 *  Earlier uploaded: https://github.com/kiselev2013/MTCalc_with_DFP_COCR/blob/master/Code/OutputNu/OutputArbitrary.h         
 *  Novosibirsk State Technical University,                                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                               
 * 
*/ 

#pragma once
#include <functional>
#include "T_Mapping.h"
#define OUTPUT
//------------------------------------------------------------------------      
// receiver with coordinates and receiver number - for sorting                  
//------------------------------------------------------------------------      
class PointRes2
{
public:
	long num;
	double point[3];
	~PointRes2() {}
};
//------------------------------------------------------------------------        
// receiver sorting predicates                                                    
//------------------------------------------------------------------------        
struct PointRes_less_x : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[0] < b.point[0];}
};
struct PointRes_less_y : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[1] < b.point[1];}
};
struct PointRes_less_z : public binary_function<PointRes2, PointRes2, bool>
{
	bool operator() (const PointRes2& a, const PointRes2& b) const {return a.point[2] < b.point[2];}
};
class Long_double_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return fabs(b.d - a.d) < 1e-6;}
};
class Long_double_num_eq : public unary_function<long_double, bool>
{
	long_double a;
public:
	explicit Long_double_num_eq(const long_double& aa) : a(aa) {}
	bool operator() (const long_double& b) const {return a.i == b.i;}
};
struct Long_double_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.d < b.d;}
};
struct Long_double_num_less : public binary_function<long_double, long_double, bool>
{
	bool operator() (const long_double& a, const long_double& b) const {return a.i < b.i;}
};
//------------------------------------------------------------------------      
// Base class for output                                                        
//------------------------------------------------------------------------      
struct Output3dArbitrary
{
	int withSpline3d;
	int withSpline2d;
	int zeroPlane;   //output on the plane z=0 through the Gauss points on the upper faces of the elements 
	long kpar;
	long kuzlov;
	long n_pointres;
	double (*pointres)[3];
	double (*xyz)[3];

	long (*nver)[14];
	long (*nvtr)[8];
	long *type_of_hex;  // 0..30 - parallelepiped, 31..61 - hexahedron 

	vector<long> elemForPoint;               // for each receiver stores the element where it will fall                    
	vector< vector<long> > PointresForElem;  // for each element stores the numbers of receivers that fall into it         
	vector<long_double> PointresXsorted;
	vector<long_double> PointresYsorted;
	vector<long_double> PointresZsorted;

	Output3dArbitrary(int withSpline3d, int withSpline2d, int zeroPlane, long kuzlov, long kpar, long n_pointres,
		double (*pointres)[3], double (*xyz)[3]);
	~Output3dArbitrary();

	void PointresSort();
	int FindElemForReceivers();

	long GetGlobalVertex(long nElem, long nLocVertex);
	long GetTypeOfElement(long nElem);
	int IsPointInsideBrick(double xm, double ym, double zm, double *x0, double *x1);
};
//------------------------------------------------------------------------    
// output from vector elements                                                
//------------------------------------------------------------------------    
struct OutputVect3d : public Output3dArbitrary
{
	long n_edges;
	long (*ed)[25];
	// constructor for vector task
	OutputVect3d(int withSpline3d, int withSpline2d, int zeroPlane,
		long kuzlov, long kpar, long n_edges, long n_pointres,
		double (*pointres)[3], double (*xyz)[3], long (*nver)[14], long (*ed)[25]);
	~OutputVect3d();
	int OutputFieldAtReceivers(double *v3, double *result, int derive);
};
//------------------------------------------------------------------------        
// output from node elements                                                      
//------------------------------------------------------------------------        
struct OutputNode3d : public Output3dArbitrary
{
	OutputNode3d(int withSpline3d, int withSpline2d, int zeroPlane,
		long kuzlov, long kpar, long n_pointres,
		double (*pointres)[3], double (*xyz)[3], long (*nvtr)[8], long *type_of_hex);
	~OutputNode3d();
};
