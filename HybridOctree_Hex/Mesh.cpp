#include "Mesh.h"

Mesh::~Mesh(void) { }

void Mesh::Initialize(int numHexa, int numPoint, int pointNumInAnElem) {
	if (numPoint == -1)
		v = new double[numHexa * pointNumInAnElem][3];
	else {
		v = new double[numPoint][3];
		r = new double[numPoint];
	}
	e = new int* [numHexa];
	for (int i = 0; i < numHexa; i++)
		e[i] = new int[pointNumInAnElem];
}

void Mesh::WriteToVtk(const char* fileName, double BOX_LENGTH_RATIO, double START_POINT[3], int elemType, bool writeNormal)
{
	FILE* dataFile = fopen(fileName, "w");
	if (NULL == dataFile)
	{
		std::cerr << "ErrorCode 0: Wrong file name" << fileName << std::endl;
		return;
	}
	fprintf(dataFile, "# vtk DataFile Version 2.0\n");
	fprintf(dataFile, "Hybrid Octree Mesh\n");
	fprintf(dataFile, "ASCII\n");
	fprintf(dataFile, "DATASET UNSTRUCTURED_GRID\n");
	fprintf(dataFile, "POINTS %i double\n", vNum);

	for (int i = 0; i < vNum; i++)
		fprintf(dataFile, "%lf %lf %lf\n", START_POINT[0] + BOX_LENGTH_RATIO * v[i][0], START_POINT[1] + BOX_LENGTH_RATIO * v[i][1], START_POINT[2] + BOX_LENGTH_RATIO * v[i][2]);

	int numInLine;
	if (elemType == 12) numInLine = 9;
	else if (elemType == 5) numInLine = 4;

	fprintf(dataFile, "CELLS %i %i\n", eNum, eNum * numInLine);

	if (elemType == 12)
		for (int i = 0; i < eNum; i++)
			fprintf(dataFile, "8 %i %i %i %i %i %i %i %i\n", e[i][0], e[i][1], e[i][2], e[i][3], e[i][4], e[i][5], e[i][6], e[i][7]);
	else if (elemType == 5)
		for (int i = 0; i < eNum; i++)
			fprintf(dataFile, "3 %i %i %i\n", e[i][0], e[i][1], e[i][2]);

	fprintf(dataFile, "CELL_TYPES %i\n", eNum);

	for (int i = 0; i < eNum; i++)
		fprintf(dataFile, "%i\n", elemType);

	if (writeNormal) {
		fprintf(dataFile, "POINT_DATA %i\n", vNum);
		fprintf(dataFile, "SCALARS curvature double 1\n");
		fprintf(dataFile, "LOOKUP_TABLE curvature\n");
		for (int i = 0; i < vNum; i++)
			fprintf(dataFile, "%lf\n", r[i]);
	}

	fclose(dataFile);
}

void Mesh::WriteToRaw(const char* fileName, double BOX_LENGTH_RATIO, double START_POINT[3])
{
	FILE* dataFile = fopen(fileName, "w");
	if (NULL == dataFile)
	{
		std::cerr << "ErrorCode 0: Wrong file name" << fileName << std::endl;
		return;
	}
	fprintf(dataFile, "%i %i\n", vNum, eNum);

	for (int i = 0; i < vNum; i++)
		fprintf(dataFile, "%lf %lf %lf 1\n", START_POINT[0] + BOX_LENGTH_RATIO * v[i][0], START_POINT[1] + BOX_LENGTH_RATIO * v[i][1], START_POINT[2] + BOX_LENGTH_RATIO * v[i][2]);

	for (int i = 0; i < eNum; i++)
		fprintf(dataFile, "%i %i %i %i %i %i %i %i\n", e[i][0], e[i][1], e[i][2], e[i][3], e[i][4], e[i][5], e[i][6], e[i][7]);
	fclose(dataFile);
}