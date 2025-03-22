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
 *  This file contains code of functions for binary I/O
 *  
 * 
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024 
*/

#include "stdafx.h"
#include "iobinary.h"

ifstream& operator>(ifstream& file, short &id)
{
	file.read((char*)&id, sizeof(short));
	return file;
}

ofstream& operator<(ofstream& file, const short &id)
{
	file.write((char*)&id, sizeof(short));
	return file;
}

ifstream& operator>(ifstream& file, int &id)
{
	file.read((char*)&id, sizeof(int));
	return file;
}

ofstream& operator<(ofstream& file, const int &id)
{
	file.write((char*)&id, sizeof(int));
	return file;
}

ifstream& operator>(ifstream& file, long &id)
{
	file.read((char*)&id, sizeof(long));
	return file;
}

ofstream& operator<(ofstream& file, const long &id)
{
	file.write((char*)&id, sizeof(long));
	return file;
}

ifstream& operator>(ifstream& file, double &id)
{
	file.read((char*)&id, sizeof(double));
	return file;
}

ofstream& operator<(ofstream& file, const double &id)
{
	file.write((char*)&id, sizeof(double));
	return file;
}

ifstream& operator>(ifstream& file, float &id)
{
	file.read((char*)&id, sizeof(float));
	return file;
}

ofstream& operator<(ofstream& file, const float &id)
{
	file.write((char*)&id, sizeof(float));
	return file;
}
