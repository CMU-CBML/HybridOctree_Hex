#include <ctime>
#include "HexGen.h"

using namespace std;

int main()
{
	const char* volumeFileName = "model.raw";
	const char* outputVolumeFileName = "modifiedTri.vtk";
	const char* OctreeFileName = "octree.vtk";
	const char* DualFullHexFileName = "dualFullHex.vtk";
	const char* DualHexFileName = "dualHex.vtk";
	const char* ProjHexFileName = "projHex.vtk";
	
	// 0/1/2/3/4
	int progress = 0;
	// 1/true: read existing file; 0/false: create new file
	bool octreeExist = (progress > 0);// want this->progress = 0;
	bool dualFullHexExist = (progress > 1);// want this->progress = 1;
	bool dualHexExist = (progress > 2);// want this->progress = 2;
	bool projHexExist = (progress > 3);// want this->progress = 3;
	
	clock_t start, finish;
	double duration;

	hexGen hexgen(VOXEL_SIZE);

	start = clock();
	
	hexgen.InitializeOctree(volumeFileName, outputVolumeFileName);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for reading surface mesh: " << duration << endl;
	start = clock();

	if (!octreeExist) {
		hexgen.ConstructOctree();
		hexgen.OutputOctree(OctreeFileName);
	}
	else
		hexgen.ReadOctree(OctreeFileName);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for constructing octree: " << duration << endl;
	start = clock();

	if (!dualFullHexExist)
		hexgen.DualFullHexMeshExtraction(DualFullHexFileName);
	else
		hexgen.ReadDualFullHex(DualFullHexFileName);
	
	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for generating dual mesh: " << duration << endl;
	start = clock();

	if (!dualHexExist)
		hexgen.RemoveOutsideElement(DualHexFileName);
	else
		hexgen.ReadDualHex(DualHexFileName);

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for extracting interior dual mesh: " << duration << endl;
	start = clock();

	if (!projHexExist)
		hexgen.ProjectToIsoSurface(ProjHexFileName);
	else
		hexgen.ReadDualHex(ProjHexFileName);// use the same function as above to store hex info to octreeMesh

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for projecting to input surface: " << duration << endl;
	cout << "Mesh generation finished";
	return 0;
}