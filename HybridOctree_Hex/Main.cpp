#include <ctime>
#include <cstdlib>
#include <limits>
#include <string>
#include "HexGen.h"

using namespace std;

int main(int argc, char** argv)
{
	int projectionIterations = 0;
	for (int arg = 1; arg < argc; ++arg) {
		const string option(argv[arg]);
		if (option == "--help") {
			cout << "Usage: HexGen [--seed N] [--iterations N]\n"
				<< "Legacy behavior is the default. N limits actual optimization steps;\n"
				<< "iterationMesh.vtk records step N. 0 is unlimited; no convergence stop.\n";
			return 0;
		}
		if ((option != "--seed" && option != "--iterations") || arg + 1 == argc) {
			cerr << "Unknown option or missing value: " << option << endl;
			return 1;
		}
		const string value(argv[++arg]);
		if (value.empty() || value.find_first_not_of("0123456789") != string::npos) {
			cerr << "Expected a nonnegative integer: " << value << endl;
			return 1;
		}
		unsigned long long number;
		try { number = stoull(value); }
		catch (const std::exception&) { cerr << "Integer out of range" << endl; return 1; }
		if (number > static_cast<unsigned long long>(numeric_limits<int>::max())) {
			cerr << "Integer out of range" << endl; return 1;
		}
		if (option == "--seed") srand(static_cast<unsigned int>(number));
		else projectionIterations = static_cast<int>(number);
	}
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
		hexgen.ProjectToIsoSurface(ProjHexFileName, projectionIterations);
	else
		hexgen.ReadDualHex(ProjHexFileName);// use the same function as above to store hex info to octreeMesh

	finish = clock();
	duration = (double)(finish - start) / CLOCKS_PER_SEC;
	cout << "Time elapsed for projecting to input surface: " << duration << endl;
	cout << (projectionIterations > 0 ? "Requested legacy optimization steps completed; output: iterationMesh.vtk" : "Mesh generation finished") << endl;
	return 0;
}
