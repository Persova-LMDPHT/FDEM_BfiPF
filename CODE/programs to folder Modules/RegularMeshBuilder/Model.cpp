/**
 *  This code is freely available under the following conditions:
 *
 *  1) The code is to be used only for non-commercial purposes.
 *  2) No changes and modifications to the code without prior permission of the developer.
 *  3) No forwarding the code to a third party without prior permission of the developer.
 *
 *                   FDEM_BfiPF 
 * The file contains basic utility functions and simple operations
 *
 *  Written by Prof. Marina G. Persova, D.Sc. Denis V. Vagin
 *  Novosibirsk State Technical University,                                                                                                                                         
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                                               
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                                                
 *  Version 2.0 December 10, 2024                                                                                                                                                    
*/

#include "Model.h"

namespace Model
{
	void CreateAirMaterial(Material &airMaterial)
	{
		airMaterial.SigmaH.constantValue = 1e-8;
		airMaterial.SigmaV.constantValue = 1e-8;
		airMaterial.SigmaNH.constantValue = 1e-8;
		airMaterial.SigmaNV.constantValue = 1e-8;
	}
}