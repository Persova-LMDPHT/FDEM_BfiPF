/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *                    
 * The file contains basic utility functions and simple operations
 *
 *  Earlier uploaded: https://github.com/kiselev2013/A-E_Formulations_3D_TimeDomain_EM_HorizontalGroundedWireSource/blob/master/CODE/programs%20to%20folder%20Modules/RegularMeshBuilder/BasicOperations.h    
 *  Novosibirsk State Technical University,                                                                                                                                                        
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                                                              
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                                                               
 *  Version 2.0 16 January, 2023                                                                                                                                                                   
*/

#pragma once

#include <vector>

const double _PI_ = 3.14159265358979;

namespace BasicOperations
{
	using namespace std;

	// Find index of interval containing value in v
	int BinarySearchInVectorSorted(vector<double> &v, double value);

	// Sort coordinates
	int InsertCoordinateIntoSortedArray(vector<double> &arr, double value, double eps);

	// Round double value with decrement
	double RoundDouble_(double value, double decrement);
}