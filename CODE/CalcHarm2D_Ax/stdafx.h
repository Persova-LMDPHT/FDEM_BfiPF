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
 *  This file contains list of include headers and global constants and functions
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin                                         
 *  Novosibirsk State Technical University,                                 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                       
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)        
 *  Version 2.0 December 10, 2024                                           
*/

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <tchar.h>
#include <direct.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <complex>
#include <algorithm>
#include <stdexcept> 
#include <iomanip>

using namespace std;

#define PI 3.14159265358979323846

const double MinRGlob = 1e-3;

#include "mkl.h"

extern ofstream logfile;

#define RETCODE_OK				0x0000
#define RETCODE_NOMEM			0x0001
#define RETCODE_NOFILE			0x0002
#define RETCODE_OUTOFRANGE		0x0004
#define RETCODE_SQFROMNEG		0x0008
#define RETCODE_DEVBYZERO		0x0010
#define RETCODE_NOTINIT			0x0020
#define RETCODE_BADFILE			0x0040
#define RETCODE_ERROR			0x0080
#define RETCODE_NOANOMALOBJECTS	0x0100
