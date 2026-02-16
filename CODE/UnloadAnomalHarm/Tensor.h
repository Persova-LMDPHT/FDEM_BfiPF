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
 *  This file contains structure for tensor electrical conductivity
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024  
*/

#pragma  once

class Tensor
{
public:
	double val[3][3];

	Tensor();
	Tensor(double d);
	bool Equal(Tensor &t, double d);
	bool NotEqual(Tensor &t, double d);
	inline void Clear();
};
