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
 *  This file contains headers for functions for binary I/O
 * 
 * 
 *  Written by D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia
 *  vdv_wk@mail.ru
 *  Version 2.0 December 10, 2024 
*/

#pragma once
#include <fstream>

using namespace std;

ifstream& operator>(ifstream& file, short &id);
ofstream& operator<(ofstream& file, const short &id);
ifstream& operator>(ifstream& file, int &id);
ofstream& operator<(ofstream& file, const int &id);
ifstream& operator>(ifstream& file, long &id);
ofstream& operator<(ofstream& file, const long &id);
ifstream& operator>(ifstream& file, double &id);
ofstream& operator<(ofstream& file, const double &id);
ifstream& operator>(ifstream& file, float &id);
ofstream& operator<(ofstream& file, const float &id);
