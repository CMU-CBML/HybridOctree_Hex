#ifndef MESH_H
#define MESH_H

#include <iostream>

class Mesh
{
public:
	double(*v)[3];// vertices
	int** e;// elements, each points to the address of a vertice
	int vNum = 0, eNum = 0;
	double(*r);// radius
	//destructor
	~Mesh();

	// distribute initial space to the mesh
	// pointNumInAnElem:
	// 8: Hexahedral (default)
	// 3: Triangular
	void Initialize(int numHexa, int numPoint = -1, int pointNumInAnElem = 8);

	// write to vtk file
	void WriteToVtk(const char* fileName, double BOX_LENGTH_RATIO, double START_POINT[3], int elemType = 12, bool writeNormal = false);

	// write to raw file
	void WriteToRaw(const char* fileName, double BOX_LENGTH_RATIO, double START_POINT[3]);
};

#endif MESH_H