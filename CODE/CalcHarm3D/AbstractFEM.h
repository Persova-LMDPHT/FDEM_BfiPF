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
 *  This file contains headers of abstract class to smooth output in 3D tasks
 *  
 *  
 *  Written by D.Sc. Denis V. Vagin, Prof. Marina G. Persova
 *  Novosibirsk State Technical University, 
 *  20 Prospekt K. Marksa, Novosibirsk,630073, Russia 
 *  Corresponding author: mpersova@mail.ru (Prof. Marina G. Persova) 
 * Version 2.0 December 10, 2024                                  
 * 
*/

#pragma once
#include "PointVector.h"

enum Res3DValueType
{
	vtWithDiscontinuity,
	vtWithoutDiscontinuity,
	vtU,
	vtdUdx,
	vtdUdy,
	vtdUdz,
	vtExSin,
	vtExCos,
	vtEySin,
	vtEyCos,
	vtEzSin,
	vtEzCos,
	vtAx,
	vtAy,
	vtAz,
	vtAxSin,
	vtAxCos,
	vtAySin,
	vtAyCos,
	vtAzSin,
	vtAzCos,
	vtRotxA,
	vtRotyA,
	vtRotzA,
	vtRotxASin,
	vtRotxACos,
	vtRotyASin,
	vtRotyACos,
	vtRotzASin,
	vtRotzACos
};

class AbstractFEM3D
{
public:
	virtual int GetNumberOfNodes()=0;
	virtual int GetNumberOfElements()=0;
	virtual int GetElementNodesNumber()=0;
	virtual const pv::Point3D GetNode(const int& i_node)=0;
	virtual const pv::Point3D GetNodeTrue(const int& i_node)=0;
	virtual int GetNodeNumberOnElement(const int& i_element, const int& i_node)=0;
	virtual int GetElementMaterial(const int& i_element)=0;
	virtual int GetTypeOfElement(const int& i_element)=0;
	virtual double GetValueInElementCenter(const int& i_element, const Res3DValueType& r_type)=0;
	virtual int GetNumberOfResPoints(const Res3DValueType& r_type)=0;
	virtual pv::Point3D GetResPoint(const Res3DValueType& r_type, const int& i_point)=0;
	virtual int * GetPointerToRegular()=0;
	virtual int GetXSize()=0;
	virtual int GetYSize()=0;
	virtual int GetZSize()=0;
	virtual double *GetPointerToX()=0;
	virtual double *GetPointerToY()=0;
	virtual double *GetPointerToZ()=0;
	virtual void SaveResult(const Res3DValueType& r_type, const double& r_value, const int& i_point, const int& i_time, int ipls)=0;
};
