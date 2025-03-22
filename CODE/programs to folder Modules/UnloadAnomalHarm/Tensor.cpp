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
 *  This file contains code for tensor electrical conductivity
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024  
*/

#include "stdafx.h"
#include "Tensor.h"

Tensor::Tensor()
{
	Clear();
}

void Tensor::Clear()
{
	val[0][0] = 0.0;
	val[0][1] = 0.0;
	val[0][2] = 0.0;

	val[1][0] = 0.0;
	val[1][1] = 0.0;
	val[1][2] = 0.0;

	val[2][0] = 0.0;
	val[2][1] = 0.0;
	val[2][2] = 0.0;
}

Tensor::Tensor(double d)
{
	val[0][0] = d;
	val[0][1] = 0.0;
	val[0][2] = 0.0;

	val[1][0] = 0.0;
	val[1][1] = d;
	val[1][2] = 0.0;

	val[2][0] = 0.0;
	val[2][1] = 0.0;
	val[2][2] = d;
}

bool Tensor::Equal(Tensor &t, double d)
{
	int i, j;
	const double eps = 1e-9;

	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++)
		{
			if (i == j)
			{
				if (fabs(t.val[i][i] - d) > eps)
					return false;
			}
			else
			{
				if (fabs(t.val[i][j]) > eps)
					return false;
			}
		}
	}

	return true;
}

bool Tensor::NotEqual(Tensor &t, double d)
{
	return !Equal(t, d);
}
