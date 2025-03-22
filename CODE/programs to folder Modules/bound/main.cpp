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
 *  This file contains the code for output FEM information about layers boundaries.
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024 
*/

#include "stdafx.h"
#include "task2d.h"

ofstream logfile;

struct layer
{
	int k;
	double sig,eps;
};

int main()
{
	int i,j,k,retp,nelmJ,kr,kz,ilay,nlay;
	ifstream inf;
	ofstream ofp;
	vector<int> ElementFlag;
	vector<layer> layers;

	logfile.open("LogOutput2dUBound");

	task2d task("Ax");
	retp=task.Read();
	if(retp)
	{
		logfile<<"Function Read() for Ax returned "<<retp<<'\n';
		return 1;
	}

	ElementFlag.resize(task.krect);
	kr=task.qr-1;
	kz=task.qz-1;

	ofp.open("ElemFlag",ios::binary);
	for(i=0;i<task.krect;i++)
	{
		ofp<0;
	}
	ofp.close();
	ofp.clear();


	for(i=0;i<task.krect;i++)ElementFlag[i]=1;

	nlay=0;
	layers.clear();

	i=task.reg[0]-1;
	for(j=1;j<kz;j++)
	{
		k=task.reg[j*kr]-1;
		if(ElementFlag[k])
		{
			if(k!=i && task.rect[k].mtr!=task.rect[i].mtr)
			{
				if(nlay)
				{
					layers.resize(nlay+1);
					layers[nlay].k=j;
					layers[nlay].sig=task.sigma[task.rect[k].mtr-1];
					layers[nlay].eps=task.dpr[task.rect[k].mtr-1];
					nlay++;
				}
				else
				{
					layers.resize(2);
					layers[0].k=0;
					layers[0].sig=task.sigma[task.rect[i].mtr-1];
					layers[0].eps=task.dpr[task.rect[i].mtr-1];
					layers[1].k=j;
					layers[1].sig=task.sigma[task.rect[k].mtr-1];
					layers[1].eps=task.dpr[task.rect[k].mtr-1];
					nlay=2;
				}
			}
			ElementFlag[k]=0;
		}
		i=k;
	}

	for(i=0;i<task.krect;i++)ElementFlag[i]=1;


	nelmJ = 0;
	ofp.open("jump.dat", ios::binary);
	for(ilay=1;ilay<nlay;ilay++)
	{
		j=layers[ilay].k;
		for(i=0;i<kr;i++)
		{
			k=task.reg[(j-1)*kr+i]-1;
			if(ElementFlag[k])
			{
				ofp<k+1<1<layers[ilay-1].sig<layers[ilay-1].eps<layers[ilay].sig<layers[ilay].eps;
				ElementFlag[k]=0;
				nelmJ++;
			}
			k=task.reg[j*kr+i]-1;
			if(ElementFlag[k])
			{
				ofp<k+1<2<layers[ilay-1].sig<layers[ilay-1].eps<layers[ilay].sig<layers[ilay].eps;
				ElementFlag[k]=0;
				nelmJ++;
			}
		}
	}
	ofp.close();
	ofp.clear();


	ofp.open("j_elm.txt");
	ofp<<nelmJ<<'\n';
	ofp.close();
	ofp.clear();


	ElementFlag.clear();
	layers.clear();

	logfile.close();
	logfile.clear();

	return 0;
}
