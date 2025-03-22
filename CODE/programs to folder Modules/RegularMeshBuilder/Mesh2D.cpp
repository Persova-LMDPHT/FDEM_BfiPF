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
 *  Earlier uploaded: https://github.com/kiselev2013/A-E_Formulations_3D_TimeDomain_EM_HorizontalGroundedWireSource/blob/master/CODE/programs%20to%20folder%20Modules/RegularMeshBuilder/Mesh2D.cpp   
 *  Novosibirsk State Technical University,                                                                                                                                              
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia                                                                                                                                    
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova)                                                                                                                     
 *  Version 2.0 16 January, 2023                                                                                                                                                         
*/

#include "Mesh2D.h"
#include "Mesh1D.h"


namespace Mesh2D
{
	int BuildBaseMesh(Settings &meshSettings, double **modelBounds, vector<double> *mesh)
	{
		try
		{
			int direction;
			direction = 0; if (Mesh1D::BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh2D : Could not build X base mesh\n"); return 1; }
			direction = 1; if (Mesh1D::BuildBaseMesh1D(meshSettings.steps[direction], meshSettings.sparse[direction], modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : BuildBaseMesh2D : Could not build Y base mesh\n"); return 1; }
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : BuildBaseMesh3D : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	int AlignMeshBounds(double **modelBounds, vector<double> *mesh)
	{
		try
		{
			int direction;
			direction = 0; if (Mesh1D::AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align X mesh\n"); return 1; }
			direction = 1; if (Mesh1D::AlignMeshBounds(modelBounds[direction], mesh[direction]) != 0) { write_to_log("Error : AlignMeshBounds : Could not align Y mesh\n"); return 1; }
			return 0;
		}
		catch (exception ex)
		{
			write_to_log("Error : AlignMeshBounds : ");
			write_to_log(ex.what());
			write_to_log("\n");
			return 1;
		}
	}

	vector<int> GetFirstBoundaryConditionNodesNumbers(int nx, int ny, int index)
	{
		vector<int> nodes;
		nodes.reserve(nx + ny);
		bool left = index == 1 || index == 2;
		bool top = index != 4;

		if (left)
			for (int nodeIndex = nx; nodeIndex < nx * ny - nx; nodeIndex += nx)
				nodes.push_back(nodeIndex + 1);

		//bottom
		for (int nodeIndex = 0; nodeIndex < nx; ++nodeIndex)
			nodes.push_back(nodeIndex + 1);

		//right
		for (int nodeIndex = 2 * nx - 1; nodeIndex < nx * ny - 1; nodeIndex += nx)
			nodes.push_back(nodeIndex);

		if (top)
			for (int nodeIndex = nx * ny - nx; nodeIndex < nx * ny; ++nodeIndex)
				nodes.push_back(nodeIndex);

		return nodes;

	}

}