// Replicates HybridOctree_Hex's two octree-refinement criteria (per-vertex
// dihedral "curvature" r[] from ReadRawData, and the per-triangle ray-cast
// "thickness" from GetCellValue) WITHOUT running the expensive
// ComputeCellValue() octree sweep, so C_THRES / H_THRES can be swept cheaply.
//
// All numerical code below is copied verbatim from HybridOctree_Hex_v1.0's
// HexGen.cpp / Initialization.h (Intersect(), the r[] dihedral accumulation in
// ReadRawData(), and the thickness loop in GetCellValue()), so the counts it
// reports are exactly the refineTri*/refineTriPt* list sizes HexGen itself
// would build for the given thresholds.
//
// Usage: refine_criteria_stats <model.raw> [c0 c1 c2 c3 c4] [h0 h1 h2 h3 h4]
// Build: c++ -std=c++17 -O2 -o refine_criteria_stats refine_criteria_stats.cpp

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <array>

static const double PI = 3.1415926535897932384626433;
static const double DIST_THRES = 1e-12;

#define CROSS(dest, v1, v2)   dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
                              dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
                              dest[2] = v1[0] * v2[1] - v1[1] * v2[0];
#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

static inline double dist(double x, double y, double z, double v2[3]) {
    return (x - v2[0]) * (x - v2[0]) + (y - v2[1]) * (y - v2[1]) + (z - v2[2]) * (z - v2[2]);
}

// verbatim from HexGen.cpp:164 (ray/triangle)
static inline int Intersect(double a[3], double b[3], double c[3], double p[3], double dir[3], double* e, double& alpha) {
    double A = (c[1] * b[2] - b[1] * c[2] + a[1] * c[2] - a[2] * c[1] - a[1] * b[2] + a[2] * b[1]);
    double B = (a[0] * (b[2] - c[2]) - b[0] * (a[2] - c[2]) + c[0] * (a[2] - b[2]));
    double C = (a[0] * (c[1] - b[1]) - b[0] * (c[1] - a[1]) + c[0] * (b[1] - a[1]));
    double D = a[0] * (b[1] * c[2] - c[1] * b[2]) - b[0] * (a[1] * c[2] - a[2] * c[1]) + c[0] * (a[1] * b[2] - a[2] * b[1]);
    double tmp2 = A * dir[0] + B * dir[1] + C * dir[2];
    if (std::abs(tmp2) < DIST_THRES) return -1;
    alpha = (-A * p[0] - B * p[1] - C * p[2] - D) / tmp2;
    e[0] = p[0] + dir[0] * alpha;
    e[1] = p[1] + dir[1] * alpha;
    e[2] = p[2] + dir[2] * alpha;
    double AP[3], AC[3], AB[3];
    AP[0] = e[0] - a[0]; AP[1] = e[1] - a[1]; AP[2] = e[2] - a[2];
    AC[0] = c[0] - a[0]; AC[1] = c[1] - a[1]; AC[2] = c[2] - a[2];
    AB[0] = b[0] - a[0]; AB[1] = b[1] - a[1]; AB[2] = b[2] - a[2];
    double fI = DOT(AP, AC) * DOT(AB, AB) - DOT(AP, AB) * DOT(AC, AB);
    double fJ = DOT(AP, AB) * DOT(AC, AC) - DOT(AP, AC) * DOT(AB, AC);
    double fD = DOT(AC, AC) * DOT(AB, AB) - DOT(AC, AB) * DOT(AC, AB);
    if (fI > 0 && fJ > 0 && fI + fJ < fD) return 1;
    if (fI == 0 || fJ == 0 || fI + fJ == fD) return -1;
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <model.raw> [c0..c4] [h0..h4]\n", argv[0]); return 1; }
    double C_THRES[5] = { 0.15, 0.3, 0.6, 1.2, 2.4 };
    double H_THRES[5] = { 16, 8, 4, 2, 1 };
    if (argc >= 7) for (int i = 0; i < 5; i++) C_THRES[i] = atof(argv[2 + i]);
    if (argc >= 12) for (int i = 0; i < 5; i++) H_THRES[i] = atof(argv[7 + i]);

    FILE* f = fopen(argv[1], "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", argv[1]); return 1; }
    char line[256];
    int points, elements;
    if (!fgets(line, sizeof(line), f) || sscanf(line, "%d %d", &points, &elements) != 2) return 1;

    std::vector<std::array<double, 3>> v(points);
    std::vector<std::array<int, 3>> e(elements);
    std::vector<double> r(points, 0.0);
    double box[3][2] = { {1e30,-1e30},{1e30,-1e30},{1e30,-1e30} };
    for (int i = 0; i < points; i++) {
        if (!fgets(line, sizeof(line), f)) { fprintf(stderr, "short file at vertex %d\n", i); return 1; }
        sscanf(line, "%lf %lf %lf", &v[i][0], &v[i][1], &v[i][2]);
        for (int j = 0; j < 3; j++) {
            box[j][0] = std::min(box[j][0], v[i][j]);
            box[j][1] = std::max(box[j][1], v[i][j]);
        }
    }
    double BOX_LENGTH = std::max({ box[0][1] - box[0][0], box[1][1] - box[1][0], box[2][1] - box[2][0] });
    double SP[3] = { 0.5 * (box[0][0] + box[0][1] - BOX_LENGTH),
                     0.5 * (box[1][0] + box[1][1] - BOX_LENGTH),
                     0.5 * (box[2][0] + box[2][1] - BOX_LENGTH) };
    for (int i = 0; i < points; i++)
        for (int j = 0; j < 3; j++) v[i][j] = (v[i][j] - SP[j]) * 100.0 / BOX_LENGTH;
    // NOTE: replicates ReadRawData's unchecked fgets — a header overstating the
    // triangle count silently duplicates the previous line, exactly as HexGen does.
    for (int i = 0; i < elements; i++) {
        fgets(line, sizeof(line), f);
        sscanf(line, "%d %d %d", &e[i][0], &e[i][1], &e[i][2]);
    }
    fclose(f);

    // curvature r[] — verbatim structure from ReadRawData
    int publicIdx[2];
    double l1[3], lPublic[3], l2[3], cross1[3], cross2[3], angle;
    // build incidence to keep it tractable (same result, faster than the O(V*E^2) original)
    std::vector<std::vector<int>> inc(points);
    for (int i = 0; i < elements; i++) for (int k = 0; k < 3; k++) inc[e[i][k]].push_back(i);
    for (int i = 0; i < points; i++) {
        r[i] = 0;
        const std::vector<int>& L = inc[i];
        for (size_t jj = 0; jj + 1 < L.size(); jj++) {
            int j = L[jj];
            int k = (e[j][0] == i) ? 0 : (e[j][1] == i ? 1 : 2);
            for (size_t ll = jj + 1; ll < L.size(); ll++) {
                int l = L[ll];
                int m = (e[l][0] == i) ? 0 : (e[l][1] == i ? 1 : 2);
                if (e[j][(k + 1) % 3] == e[l][(m + 1) % 3]) { publicIdx[0] = 1; publicIdx[1] = 1; }
                else if (e[j][(k + 1) % 3] == e[l][(m + 2) % 3]) { publicIdx[0] = 1; publicIdx[1] = 2; }
                else if (e[j][(k + 2) % 3] == e[l][(m + 1) % 3]) { publicIdx[0] = 2; publicIdx[1] = 1; }
                else if (e[j][(k + 2) % 3] == e[l][(m + 2) % 3]) { publicIdx[0] = 2; publicIdx[1] = 2; }
                else continue;
                for (int d = 0; d < 3; d++) {
                    l1[d] = -v[i][d] + v[e[j][(k + 3 - publicIdx[0]) % 3]][d];
                    l2[d] = -v[i][d] + v[e[l][(m + 3 - publicIdx[1]) % 3]][d];
                    lPublic[d] = -v[i][d] + v[e[j][(k + publicIdx[0]) % 3]][d];
                }
                CROSS(cross1, l1, lPublic)
                CROSS(cross2, l2, lPublic)
                angle = DOT(cross1, cross2) / sqrt(dist(0, 0, 0, cross1) * dist(0, 0, 0, cross2));
                angle = (angle >= -1) ? angle : -1;
                angle = (angle <= 1) ? acos(angle) : 0;
                r[i] += (angle - PI) * (angle - PI);
                break;
            }
        }
    }

    // curvature candidate lists (GetCellValue's first half)
    std::unordered_set<int> cPt[5], cTri[5];
    for (int i = 0; i < elements; i++)
        for (int j = 0; j < 3; j++)
            for (int t = 0; t < 5; t++) {
                if (r[e[i][j]] > C_THRES[t]) { cTri[t].insert(i); cPt[t].insert(e[i][j]); }
                else break;
            }

    // thickness candidate lists (GetCellValue's second half) — verbatim
    std::unordered_set<int> hPt[5], hTri[5];
    double center[3], tmp[3], dir[3], len;
    for (int i = 0; i < elements; i++) {
        for (int d = 0; d < 3; d++) { center[d] = v[e[i][1]][d] - v[e[i][0]][d]; tmp[d] = v[e[i][2]][d] - v[e[i][0]][d]; }
        CROSS(dir, center, tmp)
        len = sqrt(dist(0, 0, 0, dir));
        dir[0] /= len; dir[1] /= len; dir[2] /= len;
        for (int d = 0; d < 3; d++) center[d] = (v[e[i][0]][d] + v[e[i][1]][d] + v[e[i][2]][d]) / 3;
        for (int j = i + 1; j < elements; j++) {
            int k = Intersect(v[e[j][0]].data(), v[e[j][1]].data(), v[e[j][2]].data(), center, dir, tmp, len);
            tmp[0] = std::abs(dir[0]) > std::abs(dir[1]) ? std::abs(dir[0]) : std::abs(dir[1]);
            tmp[0] = std::abs(dir[2]) > tmp[0] ? std::abs(dir[2]) : tmp[0];
            len = tmp[0] * std::abs(len);
            if (k == 1 && len > 0.125 * H_THRES[4])
                for (int t = 0; t < 5; t++) {
                    if (len < H_THRES[t]) {
                        hTri[t].insert(i); hTri[t].insert(j);
                        for (int d = 0; d < 3; d++) hPt[t].insert(e[i][d]);
                    }
                    else break;
                }
        }
    }

    printf("model=%s  points=%d triangles=%d\n", argv[1], points, elements);
    printf("C_THRES = {%g,%g,%g,%g,%g}\n", C_THRES[0], C_THRES[1], C_THRES[2], C_THRES[3], C_THRES[4]);
    printf("H_THRES = {%g,%g,%g,%g,%g}\n", H_THRES[0], H_THRES[1], H_THRES[2], H_THRES[3], H_THRES[4]);
    printf("%-6s %-6s %-10s %-10s %-10s %-10s %-10s %-10s\n",
        "idx", "level", "curvPt", "curvTri", "thickPt", "thickTri", "unionPt", "unionTri");
    for (int t = 0; t < 5; t++) {
        std::unordered_set<int> uPt = cPt[t], uTri = cTri[t];
        for (int x : hPt[t]) uPt.insert(x);
        for (int x : hTri[t]) uTri.insert(x);
        printf("%-6d %-6d %-10zu %-10zu %-10zu %-10zu %-10zu %-10zu\n",
            t, t + 4, cPt[t].size(), cTri[t].size(), hPt[t].size(), hTri[t].size(), uPt.size(), uTri.size());
    }
    // curvature distribution, for choosing thresholds
    std::vector<double> rs(r);
    std::sort(rs.begin(), rs.end());
    printf("curvature r[] percentiles: ");
    for (int p : {50, 75, 90, 95, 97, 98, 99}) printf(" p%d=%.4f", p, rs[(size_t)(rs.size() * p / 100)]);
    printf("  max=%.4f\n", rs.back());
    return 0;
}
