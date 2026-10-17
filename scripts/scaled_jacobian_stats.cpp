// Computes min/max scaled Jacobian for a VTK ASCII UNSTRUCTURED_GRID hex mesh
// (VTK_HEXAHEDRON, cell type 12) — the metric Table 2 of Tong et al. 2024
// reports as "Scaled Jacobian [worst; best]".
//
// The Sj() function below is copied verbatim from HybridOctree_Hex/HexGen.cpp
// (lines 283-538 as of this writing) — the exact scaled-Jacobian
// implementation this codebase itself uses internally, so this tool's output
// is by construction consistent with what HexGen already computed while
// generating the mesh, not an independent reimplementation that could
// silently diverge from it via a transcription slip. It's also the same
// standard 9-sample (body-center + 8 corners) normalized-triple-product
// formulation used in the sibling HexOpt repo's meshQuality.cpp (sJGrad),
// confirmed by inspection to share the same volume/(len0*len1*len2)
// structure, just shaped for gradient computation there instead of a pure
// value. See analysis.md's summary log (2026-08-04 entry) for that
// cross-check.
//
// Usage: scaled_jacobian_stats <mesh.vtk>
// Build: c++ -std=c++17 -O2 -o scaled_jacobian_stats scaled_jacobian_stats.cpp

#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Constants from HybridOctree_Hex/Initialization.h, needed by Sj() below.
static const int MAX_NUM2 = 2147483647;
static const double DIST_THRES = 1e-12;// threshold when judging point overlap

// --- verbatim copy of hexGen::Sj() from HexGen.cpp begins ---
inline double Sj(double p0[3], double p1[3], double p2[3], double p3[3], double p4[3], double p5[3], double p6[3], double p7[3]) {
	double minNormSJ = MAX_NUM2;

	double x0 = p1[0] + p2[0] + p5[0] + p6[0] - p0[0] - p3[0] - p4[0] - p7[0];
	double y0 = p1[1] + p2[1] + p5[1] + p6[1] - p0[1] - p3[1] - p4[1] - p7[1];
	double z0 = p1[2] + p2[2] + p5[2] + p6[2] - p0[2] - p3[2] - p4[2] - p7[2];

	double x1 = p2[0] + p3[0] + p6[0] + p7[0] - p0[0] - p1[0] - p4[0] - p5[0];
	double y1 = p2[1] + p3[1] + p6[1] + p7[1] - p0[1] - p1[1] - p4[1] - p5[1];
	double z1 = p2[2] + p3[2] + p6[2] + p7[2] - p0[2] - p1[2] - p4[2] - p5[2];

	double x2 = p4[0] + p5[0] + p6[0] + p7[0] - p0[0] - p1[0] - p2[0] - p3[0];
	double y2 = p4[1] + p5[1] + p6[1] + p7[1] - p0[1] - p1[1] - p2[1] - p3[1];
	double z2 = p4[2] + p5[2] + p6[2] + p7[2] - p0[2] - p1[2] - p2[2] - p3[2];

	double volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	double len0 = x0 * x0 + y0 * y0 + z0 * z0;
	double len1 = x1 * x1 + y1 * y1 + z1 * z1;
	double len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	double len = std::sqrt(len0 * len1 * len2);
	double tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(0,0,0)
	x0 = p1[0] - p0[0];
	y0 = p1[1] - p0[1];
	z0 = p1[2] - p0[2];

	x1 = p3[0] - p0[0];
	y1 = p3[1] - p0[1];
	z1 = p3[2] - p0[2];

	x2 = p4[0] - p0[0];
	y2 = p4[1] - p0[1];
	z2 = p4[2] - p0[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(1,0,0)
	x0 = p2[0] - p1[0];
	y0 = p2[1] - p1[1];
	z0 = p2[2] - p1[2];

	x1 = p0[0] - p1[0];
	y1 = p0[1] - p1[1];
	z1 = p0[2] - p1[2];

	x2 = p5[0] - p1[0];
	y2 = p5[1] - p1[1];
	z2 = p5[2] - p1[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(1,1,0)
	x0 = p3[0] - p2[0];
	y0 = p3[1] - p2[1];
	z0 = p3[2] - p2[2];

	x1 = p1[0] - p2[0];
	y1 = p1[1] - p2[1];
	z1 = p1[2] - p2[2];

	x2 = p6[0] - p2[0];
	y2 = p6[1] - p2[1];
	z2 = p6[2] - p2[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(0,1,0)
	x0 = p0[0] - p3[0];
	y0 = p0[1] - p3[1];
	z0 = p0[2] - p3[2];

	x1 = p2[0] - p3[0];
	y1 = p2[1] - p3[1];
	z1 = p2[2] - p3[2];

	x2 = p7[0] - p3[0];
	y2 = p7[1] - p3[1];
	z2 = p7[2] - p3[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(0,0,1)
	x0 = p7[0] - p4[0];
	y0 = p7[1] - p4[1];
	z0 = p7[2] - p4[2];

	x1 = p5[0] - p4[0];
	y1 = p5[1] - p4[1];
	z1 = p5[2] - p4[2];

	x2 = p0[0] - p4[0];
	y2 = p0[1] - p4[1];
	z2 = p0[2] - p4[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(1,0,1)
	x0 = p4[0] - p5[0];
	y0 = p4[1] - p5[1];
	z0 = p4[2] - p5[2];

	x1 = p6[0] - p5[0];
	y1 = p6[1] - p5[1];
	z1 = p6[2] - p5[2];

	x2 = p1[0] - p5[0];
	y2 = p1[1] - p5[1];
	z2 = p1[2] - p5[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(1,1,1)
	x0 = p5[0] - p6[0];
	y0 = p5[1] - p6[1];
	z0 = p5[2] - p6[2];

	x1 = p7[0] - p6[0];
	y1 = p7[1] - p6[1];
	z1 = p7[2] - p6[2];

	x2 = p2[0] - p6[0];
	y2 = p2[1] - p6[1];
	z2 = p2[2] - p6[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		minNormSJ = tempNormSJ;

	// J(0,1,1)
	x0 = p6[0] - p7[0];
	y0 = p6[1] - p7[1];
	z0 = p6[2] - p7[2];

	x1 = p4[0] - p7[0];
	y1 = p4[1] - p7[1];
	z1 = p4[2] - p7[2];

	x2 = p3[0] - p7[0];
	y2 = p3[1] - p7[1];
	z2 = p3[2] - p7[2];

	volume = x0 * (y1 * z2 - y2 * z1) + y0 * (z1 * x2 - z2 * x1) + z0 * (x1 * y2 - x2 * y1);

	len0 = x0 * x0 + y0 * y0 + z0 * z0;
	len1 = x1 * x1 + y1 * y1 + z1 * z1;
	len2 = x2 * x2 + y2 * y2 + z2 * z2;

	if (len0 <= DIST_THRES || len1 <= DIST_THRES || len2 <= DIST_THRES)
		return -(double)MAX_NUM2;

	len = std::sqrt(len0 * len1 * len2);
	tempNormSJ = volume / len;

	if (tempNormSJ < minNormSJ)
		return tempNormSJ;
	return minNormSJ;
}
// --- verbatim copy of hexGen::Sj() ends ---

int main(int argc, char** argv) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <mesh.vtk>" << std::endl;
		return 1;
	}
	std::ifstream in(argv[1]);
	if (!in) {
		std::cerr << "Cannot open " << argv[1] << std::endl;
		return 1;
	}

	std::string line;
	int pointNum = 0, cellNum = 0;
	std::vector<std::array<double, 3>> pts;

	// Skip header, find POINTS.
	while (std::getline(in, line)) {
		if (line.rfind("POINTS", 0) == 0) {
			std::istringstream iss(line);
			std::string tag, type;
			iss >> tag >> pointNum >> type;
			break;
		}
	}
	pts.resize(pointNum);
	for (int i = 0; i < pointNum; ++i)
		in >> pts[i][0] >> pts[i][1] >> pts[i][2];

	// Find CELLS.
	while (std::getline(in, line)) {
		if (line.rfind("CELLS", 0) == 0) {
			std::istringstream iss(line);
			std::string tag;
			int totalInts;
			iss >> tag >> cellNum >> totalInts;
			break;
		}
	}

	double worst = std::numeric_limits<double>::infinity();
	double best = -std::numeric_limits<double>::infinity();
	long long hexCount = 0, nonHexSkipped = 0;

	for (int c = 0; c < cellNum; ++c) {
		int n;
		in >> n;
		std::vector<int> idx(n);
		for (int i = 0; i < n; ++i) in >> idx[i];
		if (n != 8) {// only VTK_HEXAHEDRON (type 12, 8 points) contributes to
		             // Table 2's scaled-Jacobian range
			++nonHexSkipped;
			continue;
		}
		double p[8][3];
		for (int i = 0; i < 8; ++i) {
			p[i][0] = pts[idx[i]][0];
			p[i][1] = pts[idx[i]][1];
			p[i][2] = pts[idx[i]][2];
		}
		double sj = Sj(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]);
		if (sj < worst) worst = sj;
		if (sj > best) best = sj;
		++hexCount;
	}

	std::printf("Points:   %d\n", pointNum);
	std::printf("Cells:    %d (%lld hex, %lld non-hex skipped)\n", cellNum, hexCount, nonHexSkipped);
	std::printf("Worst SJ: %.6f\n", worst);
	std::printf("Best SJ:  %.6f\n", best);
	return 0;
}
