/*
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

#pragma once

#include <vector>
#include <stdio.h>
#include <fstream>

#include "Mesh3D.h"
#include "Mesh2D.h"
#include "Model.h"

namespace ReadWrite
{
	using namespace std;

	// Read int value
	int ReadInt(char *file_name, int &val);

	// Read objects from file
	int ReadObjects(char *file_name, vector<Model::GeoBox<double>> &objects);

	// Read layers from file
	int ReadZSig2D(char *file_name, vector<Model::GeoLayer<double>> &layers);

	// Read receivers from file
	int ReadReceivers(char *file_name, vector<Model::Receiver> &receivers);

	// Read generators from file
	int ReadGenerators(char *file_name, vector<Model::Generator> &generators);

	// Read settings from file
	int ReadSettings(char *file_name, Mesh3D::Settings &mesh3DSettings, Mesh2D::Settings &mesh2DSettings);

	
	// Write all mesh files procedure
	int WriteMesh3D(vector<double> *meshes1D, vector<double> *meshesTemplate1D, vector<int> &materialNumbers, set<pair<int, int>> &edges, vector<vector<int>> &elementsByEdges);

	// Write material parameters files
	int WriteMaterials3D(vector<Model::Material> &materials);

	// Write material parameters files
	int WriteMaterials2D(vector<Model::Material> &materials);

	// Write all mesh files procedure
	int WriteMesh2D(vector<double> *meshes1D, vector<int> &materialNumbers, int index);

	// Write receiver points
	int WritexyzVectorE(char *fileName, vector<Model::Receiver> &receivers, int refineStepsCount = 1);

	// Write receivers coordinates count
	int WriteRecvsE(char *fileName, vector<Model::Receiver> &receivers, int numberOfMNDipoles, int generatorsCount);

	// Write group info
	int WriteGroup(char *fileName, vector<Model::Generator> &generators);

	// Write gens
	int WriteGens(char *fileName, vector<Model::Generator> &generators);

	// Write time mesh files
	int WriteCurrentFunction(char *file_name, vector<double> &mesh);
	int WriteDeltaFunction(char *file_name, vector<double> &mesh);
	int WriteInfite0(char *file_name, vector<double> &mesh);
}