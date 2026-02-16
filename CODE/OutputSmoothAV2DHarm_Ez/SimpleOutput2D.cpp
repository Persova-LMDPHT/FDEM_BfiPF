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
 *  This file contains code for direct FEM output E and B fields in receivers.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024     
*/

#include "stdafx.h"
#include "SimpleOutput2D.h"
#include "task2d.h"

extern ofstream logfile;

// Direct output solution in receivers
void Output_Field(int nrec, double *v2es, double *v2ec, double *results, double *resultc, RectOutData *vRecOutData)
{
	double sigma;
	int k, node0, node1, node2, node3;
	double psi, eta, psi2, eta2;

	for (k = 0; k<nrec; k++)
	{
		RectOutData &rod = vRecOutData[k];

		node0 = rod.node0;
		node1 = rod.node1;
		node2 = rod.node2;
		node3 = rod.node3;

		psi = rod.psi;
		eta = rod.eta;

		psi2 = rod.psi2;
		eta2 = rod.eta2;

		sigma = rod.sigma;

		results[k] = psi2*eta2*v2es[node0] + psi*eta2*v2es[node1] + psi2*eta*v2es[node2] + psi*eta*v2es[node3];
		resultc[k] = psi2*eta2*v2ec[node0] + psi*eta2*v2ec[node1] + psi2*eta*v2ec[node2] + psi*eta*v2ec[node3];
	}
}

int FindElementsForReceivers(int nrec, int qr, int qz, double *rm, double *zm, int *reg, PointRZ *rec, RectOutData *vRecOutData, 
	PointRZ *pnt, Rect *rect, double *sigma, double rmin, double  zmin, double  rmax, double  zmax)
{
	bool fstop;
	int i, j, k, l, m;
	int imin, imax, jmin, jmax;


	imin = FindIntervalInDoubleMas(rm, qr, rmin);
	imax = FindIntervalInDoubleMas(rm, qr, rmax);
	jmin = FindIntervalInDoubleMas(zm, qz, zmin);
	jmax = FindIntervalInDoubleMas(zm, qz, zmax);

	if (imin == -1)
	{
		logfile << "No element for rmin= " << rmin << '\n';
		return 1;
	}

	if (imax == -1)
	{
		logfile << "No element for rmax= " << rmax << '\n';
		return 1;
	}

	if (jmin == -1)
	{
		logfile << "No element for zmin= " << zmin << '\n';
		return 1;
	}

	if (jmax == -1)
	{
		logfile << "No element for zmax= " << zmax << '\n';
		return 1;
	}

	fstop = false;
	for (k = 0; (k < nrec && !fstop); k++)
	{
		RectOutData &rod = vRecOutData[k];
		i = FindIntervalInDoubleMas(rm, imin, imax+1, rec[k].r);
		j = FindIntervalInDoubleMas(zm, jmin, jmax+1, rec[k].z);
		if (i != -1 && j != -1)
		{
			m = j*(qr - 1) + i;
			l = reg[m] - 1;

			rod.node0 = rect[l].nodes[0] - 1;
			rod.node1 = rect[l].nodes[1] - 1;
			rod.node2 = rect[l].nodes[2] - 1;
			rod.node3 = rect[l].nodes[3] - 1;

			rod.hr = pnt[rod.node3].r - pnt[rod.node0].r;
			rod.hz = pnt[rod.node3].z - pnt[rod.node0].z;
			rod.psi = (rec[k].r - pnt[rod.node0].r) / rod.hr;
			rod.eta = (rec[k].z - pnt[rod.node0].z) / rod.hz;
			rod.psi2 = 1.0 - rod.psi;
			rod.eta2 = 1.0 - rod.eta;

			rod.sigma = sigma[rect[l].mtr-1];
		}
		else
		{
			logfile << "No element for reciver " << k + 1 << " r= " << rec[k].r << " z= " << rec[k].z << '\n';
			fstop = true;
		}
	}
	return fstop;
}
