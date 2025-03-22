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
 *  This file contains global constants and functions
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024
*/

#include "stdafx.h"

void reverse_vec3(double *a)
{
	a[0]*=-1;
	a[1]*=-1;
	a[2]*=-1;
}

double norma3(double *a)
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

void normalize3(double *a)
{
	double s;
	s=a[0]*a[0]+a[1]*a[1]+a[2]*a[2];
	if(s)
	{
		s=sqrt(s);
		a[0]/=s;
		a[1]/=s;
		a[2]/=s;
	}
}

double mult_scal3(double *a,double *b)
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

void mult_vec3(double *a,double *b,double *c)
{
	c[0]=a[1]*b[2]-b[1]*a[2];
	c[1]=b[0]*a[2]-a[0]*b[2];
	c[2]=a[0]*b[1]-b[0]*a[1];
}

void mult_matrix33(double a[][3],double b[][3],double c[][3])
{
	int i,j,k;
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			c[i][j]=0.0;
			for(k=0;k<3;k++)
			{
				c[i][j]+=a[i][k]*b[k][j];
			}
		}
	}
}

double GetDeterminant33(double m[][3])
{
	return 
		m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1])+
		m[0][1]*(m[1][2]*m[2][0]-m[1][0]*m[2][2])+
		m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);
}

void TransposeMatrix33(double a[][3])
{
	int i,j;
	double t;
	for(i=0;i<3;i++)
	{
		for(j=0;j<i;j++)
		{
			t=a[i][j];
			a[i][j]=a[j][i];
			a[j][i]=t;
		}
	}
}

void GetRotationMatrix33(double M_T[][3],double *v,double angle)
{
	double st,ct;

	st=sin(angle);
	ct=cos(angle);
	
	M_T[0][0]=ct+(1-ct)*v[0]*v[0];
	M_T[0][1]=(1-ct)*v[0]*v[1]-st*v[2];
	M_T[0][2]=(1-ct)*v[0]*v[2]+st*v[1];

	M_T[1][0]=(1-ct)*v[1]*v[0]+st*v[2];
	M_T[1][1]=ct+(1-ct)*v[1]*v[1];
	M_T[1][2]=(1-ct)*v[1]*v[2]-st*v[0];

	M_T[2][0]=(1-ct)*v[2]*v[0]-st*v[1];
	M_T[2][1]=(1-ct)*v[2]*v[1]+st*v[0];
	M_T[2][2]=ct+(1-ct)*v[2]*v[2];
}

void InverseMatrix33(double m[][3],double m_1[][3])
{
	double det=GetDeterminant33(m);
	m_1[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;
	m_1[1][0]=-(m[1][0]*m[2][2]-m[1][2]*m[2][0])/det;
	m_1[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;
	m_1[0][1]=-(m[0][1]*m[2][2]-m[0][2]*m[2][1])/det;
	m_1[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;
	m_1[2][1]=-(m[0][0]*m[2][1]-m[0][1]*m[2][0])/det;
	m_1[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;
	m_1[1][2]=-(m[0][0]*m[1][2]-m[0][2]*m[1][0])/det;
	m_1[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;
}

void InverseMatrix33WithDet(double m[][3],double m_1[][3],double det)
{
	m_1[0][0]=(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det;
	m_1[1][0]=-(m[1][0]*m[2][2]-m[1][2]*m[2][0])/det;
	m_1[2][0]=(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det;
	m_1[0][1]=-(m[0][1]*m[2][2]-m[0][2]*m[2][1])/det;
	m_1[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])/det;
	m_1[2][1]=-(m[0][0]*m[2][1]-m[0][1]*m[2][0])/det;
	m_1[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])/det;
	m_1[1][2]=-(m[0][0]*m[1][2]-m[0][2]*m[1][0])/det;
	m_1[2][2]=(m[0][0]*m[1][1]-m[0][1]*m[1][0])/det;
}
