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
 *  This file contains headers for smooth FEM output E and B fields functions.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#pragma once
#include "task2d.h"

struct Output2D_PAR
{
	int fl;
	int lem[4];	
	double coef[4][4];
	double ves;	
	double v[2][3];
};

struct Output2D
{
	Output2D_PAR oror;
	Output2D_PAR oz;
	Output2D_PAR orfz[4],ozfr[4];
	PointRZ lrecr[4],lrecz[4];
	int llemr[4],llemz[4];
};

struct SmoothOutput2D
{
	int nrec,kpnt,krect;
	PointRZ *rec;
	PointRZ *pnt;
	Rect *rect;	
	vector<int> RecToElem;
	double *sigma;
	int nreg;
	int qr,qz;
	int *reg;
	double *rm,*zm;

	Output2D *OutRecList;

	SmoothOutput2D();
	~SmoothOutput2D();
	int Init(int p_nrec,int p_kpnt,int p_krect,PointRZ *p_rec,PointRZ *p_pnt,Rect *p_rect,int MatFlag,int *matrec,double *p_sigma,
		int p_nreg,int p_qr,int p_qz,int *p_reg,double *p_rm,double *p_zm);
	int Init(int rec_qr,int rec_qz, int p_kpnt, int p_krect, PointRZ *p_rec, PointRZ *p_pnt, Rect *p_rect, int MatFlag, int *matrec, double *p_sigma,
		int p_nreg, int p_qr, int p_qz, int *p_reg, double *p_rm, double *p_zm);
	void Output_dAzdr(double *v2s,double *v2c,double *results,double *resultc,int OutType,vector<int> &RecToSourceB);
	void Output_Ez(double *v2as,double *v2ac,double *v2vs,double *v2vc,double *results,double *resultc,int OutType,vector<int> &RecToSourceE);
	void GetSmoothForR(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int rlem[4],double rcoef[4][4],double rv[2][3],double &ves,int &fl,int stype);
	void GetSmoothForZ(int i,int j,int k,int ef,int l,int MatFlag,PointRZ &reck,int zlem[4],double zcoef[4][4],double zv[2][3],double &ves,int &fl,int stype);
	double GetSmoothValForZ(int k,int j,double *v2,int OutType);
	double GetSmoothDzForZ(int k,int j,double *v2,int OutType);

	void Output_Ez(double *v2as,double *v2ac,double *v2vs,double *v2vc,double *results,double *resultc,int OutType);
};
