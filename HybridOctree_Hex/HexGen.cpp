#include "HexGen.h"
double ELEM_THRES = 0.01;// element floating quality threshold
int octreeENum; int hexMeshENum;

// general functions
inline double dist(double v[3], double v2[3])
{
	return (v[0] - v2[0]) * (v[0] - v2[0]) + (v[1] - v2[1]) * (v[1] - v2[1]) + (v[2] - v2[2]) * (v[2] - v2[2]);
}
inline double dist(double x, double y, double z, double v2[3])
{
	return (x - v2[0]) * (x - v2[0]) + (y - v2[1]) * (y - v2[1]) + (z - v2[2]) * (z - v2[2]);
}

inline int Coplanar(double p1[3], double q1[3], double r1[3], double p2[3], double q2[3], double r2[3], double normal1[3], double normal2[3])
{
	double P1[2], Q1[2], R1[2];
	double P2[2], Q2[2], R2[2];

	double nX = ((normal1[0] < 0) ? -normal1[0] : normal1[0]),
		nY = ((normal1[1] < 0) ? -normal1[1] : normal1[1]),
		nZ = ((normal1[2] < 0) ? -normal1[2] : normal1[2]);

	// project onto plane YZ
	if ((nX > nZ) && (nX >= nY)) {
		P1[0] = q1[2]; P1[1] = q1[1];
		Q1[0] = p1[2]; Q1[1] = p1[1];
		R1[0] = r1[2]; R1[1] = r1[1];

		P2[0] = q2[2]; P2[1] = q2[1];
		Q2[0] = p2[2]; Q2[1] = p2[1];
		R2[0] = r2[2]; R2[1] = r2[1];
	}// project onto plane XZ
	else if ((nY > nZ) && (nY >= nX)) {
		P1[0] = q1[0]; P1[1] = q1[2];
		Q1[0] = p1[0]; Q1[1] = p1[2];
		R1[0] = r1[0]; R1[1] = r1[2];

		P2[0] = q2[0]; P2[1] = q2[2];
		Q2[0] = p2[0]; Q2[1] = p2[2];
		R2[0] = r2[0]; R2[1] = r2[2];
	}// project onto plane XY
	else {
		P1[0] = p1[0]; P1[1] = p1[1];
		Q1[0] = q1[0]; Q1[1] = q1[1];
		R1[0] = r1[0]; R1[1] = r1[1];

		P2[0] = p2[0]; P2[1] = p2[1];
		Q2[0] = q2[0]; Q2[1] = q2[1];
		R2[0] = r2[0]; R2[1] = r2[1];
	}
	return Intersect2D(P1, Q1, R1, P2, Q2, R2);
}
inline int CounterClockwiseIntersection2D(double p1[2], double q1[2], double r1[2], double p2[2], double q2[2], double r2[2])
{
	if (ORIENT_2D(p2, q2, p1) >= 0) {
		if (ORIENT_2D(q2, r2, p1) >= 0) {
			if (ORIENT_2D(r2, p2, p1) >= 0) return 1;
			INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
		}
		else {
			if (ORIENT_2D(r2, p2, p1) >= 0) INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
			else                            INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
		}
	}
	else {
		if (ORIENT_2D(q2, r2, p1) >= 0) {
			if (ORIENT_2D(r2, p2, p1) >= 0) INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
			else                            INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
		}
		else {
			INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
		}
	}
}
inline int Intersect2D(double p1[2], double q1[2], double r1[2], double p2[2], double q2[2], double r2[2])
{
	if (ORIENT_2D(p1, q1, r1) < 0)
		if (ORIENT_2D(p2, q2, r2) < 0)
			return CounterClockwiseIntersection2D(p1, r1, q1, p2, r2, q2);
		else
			return CounterClockwiseIntersection2D(p1, r1, q1, p2, q2, r2);
	else
		if (ORIENT_2D(p2, q2, r2) < 0)
			return CounterClockwiseIntersection2D(p1, q1, r1, p2, r2, q2);
		else
			return CounterClockwiseIntersection2D(p1, q1, r1, p2, q2, r2);
}
inline int Intersect(double p1[3], double q1[3], double r1[3], double p2[3], double q2[3], double r2[3])
{
	double dp1, dq1, dr1, dp2, dq2, dr2;
	double v1[3], v2[3];
	double N1[3], N2[3];

	// compute distance signs of p1, q1 and r1 to the plane of triangle(p2, q2, r2)
	SUB(v1, p2, r2)
		SUB(v2, q2, r2)
		CROSS(N2, v1, v2)

		SUB(v1, p1, r2)
		dp1 = DOT(v1, N2);
	SUB(v1, q1, r2)
		dq1 = DOT(v1, N2);
	SUB(v1, r1, r2)
		dr1 = DOT(v1, N2);

	if (dp1 * dq1 > 0 && dp1 * dr1 > 0)  return 0;

	// compute distance signs of p2, q2 and r2 to the plane of triangle(p1, q1, r1)
	SUB(v1, q1, p1)
		SUB(v2, r1, p1)
		CROSS(N1, v1, v2)

		SUB(v1, p2, r1)
		dp2 = DOT(v1, N1);
	SUB(v1, q2, r1)
		dq2 = DOT(v1, N1);
	SUB(v1, r2, r1)
		dr2 = DOT(v1, N1);

	if (dp2 * dq2 > 0 && dp2 * dr2 > 0) return 0;

	// permutation in a canonical form of T1's vertices
	if (abs(dp1) <= DIST_THRES) dp1 = 0;
	if (abs(dq1) <= DIST_THRES) dq1 = 0;
	if (abs(dr1) <= DIST_THRES) dr1 = 0;

	if (abs(dp2) <= DIST_THRES) dp2 = 0;
	if (abs(dq2) <= DIST_THRES) dq2 = 0;
	if (abs(dr2) <= DIST_THRES) dr2 = 0;

	if (dp1 > 0) {
		if (dq1 > 0) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
		else if (dr1 > 0) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
		else              TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
	}
	else if (dp1 < 0) {
		if (dq1 < 0) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
		else if (dr1 < 0) TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
		else              TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
	}
	else {
		if (dq1 < 0) {
			if (dr1 >= 0) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
			else          TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
		}
		else if (dq1 > 0) {
			if (dr1 > 0) TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
			else         TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
		}
		else {
			if (dr1 > 0) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
			else if (dr1 < 0) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
			else return Coplanar(p1, q1, r1, p2, q2, r2, N1, N2);
		}
	}
}

inline double TriArea(double a, double b, double c)
{
	double s = (a + b + c) / 2, area = s * (s - a) * (s - b) * (s - c);
	return sqrt((area < 0 ? 0 : area));
}
inline int Intersect(double a[3], double b[3], double c[3], double p[3], double dir[3], double* e, double& alpha) {
	double A = (c[1] * b[2] - b[1] * c[2] + a[1] * c[2] - a[2] * c[1] - a[1] * b[2] + a[2] * b[1]);
	double B = (a[0] * (b[2] - c[2]) - b[0] * (a[2] - c[2]) + c[0] * (a[2] - b[2]));
	double C = (a[0] * (c[1] - b[1]) - b[0] * (c[1] - a[1]) + c[0] * (b[1] - a[1]));
	double D = a[0] * (b[1] * c[2] - c[1] * b[2]) - b[0] * (a[1] * c[2] - a[2] * c[1]) + c[0] * (a[1] * b[2] - a[2] * b[1]);
	double tmp2 = A * dir[0] + B * dir[1] + C * dir[2];
	if (abs(tmp2) < DIST_THRES)// parallel
		return -1;
	else {
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
}
inline double PointToTri(double a[3], double b[3], double c[3], double p[3], double* q, double currMin) {
	double A = (c[1] * b[2] - b[1] * c[2] + a[1] * c[2] - a[2] * c[1] - a[1] * b[2] + a[2] * b[1]);
	double B = (a[0] * (b[2] - c[2]) - b[0] * (a[2] - c[2]) + c[0] * (a[2] - b[2]));
	double C = (a[0] * (c[1] - b[1]) - b[0] * (c[1] - a[1]) + c[0] * (b[1] - a[1]));
	double D = a[0] * (b[1] * c[2] - c[1] * b[2]) - b[0] * (a[1] * c[2] - a[2] * c[1]) + c[0] * (a[1] * b[2] - a[2] * b[1]);
	double sum = A * A + B * B + C * C, tmp = (-A * p[0] - B * p[1] - C * p[2] - D);
	double alpha = tmp / sqrt(sum);// distance 1
	if (abs(alpha) >= currMin) return abs(alpha);
	q[0] = p[0] + A * tmp / sum;
	q[1] = p[1] + B * tmp / sum;
	q[2] = p[2] + C * tmp / sum;
	alpha = abs(alpha);

	double QA = sqrt(dist(q, a)), QB = sqrt(dist(q, b)), QC = sqrt(dist(q, c)), AB = sqrt(dist(a, b)), AC = sqrt(dist(a, c)), BC = sqrt(dist(b, c));
	double S1 = TriArea(QA, QB, AB), S2 = TriArea(QA, QC, AC), S3 = TriArea(QB, QC, BC);

	double AP[3], ACl[3], ABl[3];
	AP[0] = p[0] - a[0]; AP[1] = p[1] - a[1]; AP[2] = p[2] - a[2];
	ACl[0] = c[0] - a[0]; ACl[1] = c[1] - a[1]; ACl[2] = c[2] - a[2];
	ABl[0] = b[0] - a[0]; ABl[1] = b[1] - a[1]; ABl[2] = b[2] - a[2];
	double fI = DOT(AP, ACl) * DOT(ABl, ABl) - DOT(AP, ABl) * DOT(ACl, ABl);
	double fJ = DOT(AP, ABl) * DOT(ACl, ACl) - DOT(AP, ACl) * DOT(ABl, ACl);
	double fD = DOT(ACl, ACl) * DOT(ABl, ABl) - DOT(ACl, ABl) * DOT(ACl, ABl);
	if (fI >= 0 && fJ >= 0 && fI + fJ <= fD) return alpha;// projection point inside tri
	double QtAB = S1 * 2 / AB;
	double QtBC = S3 * 2 / BC;
	double QtCA = S2 * 2 / AC;
	double beta;
	int projPointType;
	// 0: Q = A
	// 1: Q = B
	// 2: Q = C
	// 3: Q = QtAB
	// 4: Q = QtBC
	// 5: Q = QtCA
	if (QA < QB) {
		projPointType = 0;
		beta = QA;
	}
	else {
		projPointType = 1;
		beta = QB;
	}
	if (QC < beta) {
		projPointType = 2;
		beta = QC;
	}
	double k[3];// AB BC CA
	k[0] = ((b[0] - a[0]) * (q[0] - a[0]) + (b[1] - a[1]) * (q[1] - a[1]) + (b[2] - a[2]) * (q[2] - a[2])) / dist(a, b);
	k[1] = ((c[0] - b[0]) * (q[0] - b[0]) + (c[1] - b[1]) * (q[1] - b[1]) + (c[2] - b[2]) * (q[2] - b[2])) / dist(b, c);
	k[2] = ((a[0] - c[0]) * (q[0] - c[0]) + (a[1] - c[1]) * (q[1] - c[1]) + (a[2] - c[2]) * (q[2] - c[2])) / dist(c, a);
	if (k[0] > 0 && k[0] < 1 && QtAB < beta) {
		projPointType = 3;
		beta = QtAB;
	}
	if (k[1] > 0 && k[1] < 1 && QtBC < beta) {
		projPointType = 4;
		beta = QtBC;
	}
	if (k[2] > 0 && k[2] < 1 && QtCA < beta) {
		projPointType = 5;
		beta = QtCA;
	}
	switch (projPointType) {
	case 0:
		q[0] = a[0]; q[1] = a[1]; q[2] = a[2];
		break;
	case 1:
		q[0] = b[0]; q[1] = b[1]; q[2] = b[2];
		break;
	case 2:
		q[0] = c[0]; q[1] = c[1]; q[2] = c[2];
		break;
	case 3:
		q[0] = (1 - k[0]) * a[0] + k[0] * b[0];
		q[1] = (1 - k[0]) * a[1] + k[0] * b[1];
		q[2] = (1 - k[0]) * a[2] + k[0] * b[2];
		break;
	case 4:
		q[0] = (1 - k[1]) * b[0] + k[1] * c[0];
		q[1] = (1 - k[1]) * b[1] + k[1] * c[1];
		q[2] = (1 - k[1]) * b[2] + k[1] * c[2];
		break;
	case 5:
		q[0] = (1 - k[2]) * c[0] + k[2] * a[0];
		q[1] = (1 - k[2]) * c[1] + k[2] * a[1];
		q[2] = (1 - k[2]) * c[2] + k[2] * a[2];
	}
	return sqrt(alpha * alpha + beta * beta);
}

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

inline void iSj(double p0[3], double p1[3], double p2[3], double p3[3], double p4[3], double p5[3], double p6[3], double p7[3], int* minIdx) {
	minIdx[0] = -1; minIdx[1] = -1; minIdx[2] = -1; minIdx[3] = -1; minIdx[4] = -1; minIdx[5] = -1; minIdx[6] = -1; minIdx[7] = -1; minIdx[8] = -1;
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
	double len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[0] = 0;
	else {
		if (volume <= 0) minIdx[0] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[0] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[1] = 0;
	else {
		if (volume <= 0) minIdx[1] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[1] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[2] = 0;
	else {
		if (volume <= 0) minIdx[2] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[2] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[3] = 0;
	else {
		if (volume <= 0) minIdx[3] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[3] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[4] = 0;
	else {
		if (volume <= 0) minIdx[4] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[4] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[5] = 0;
	else {
		if (volume <= 0) minIdx[5] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[5] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[6] = 0;
	else {
		if (volume <= 0) minIdx[6] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[6] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[7] = 0;
	else {
		if (volume <= 0) minIdx[7] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[7] = 1;
	}

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
	len = std::sqrt(len0 * len1 * len2);

	if (len0 <= DIST_THRES3 || len1 <= DIST_THRES3 || len2 <= DIST_THRES3) minIdx[8] = 0;
	else {
		if (volume <= 0) minIdx[8] = 0;
		else if (volume <= ELEM_THRES * len) minIdx[8] = 1;
	}
}

// hexGen
hexGen::hexGen(int depth) {
	octreeDepth = depth;
	voxelSize = 1 << depth;// to the power of 2
}

hexGen::~hexGen(void) { }

void hexGen::InitializeOctree(const char* inputFileName, const char* outputFileName) {
	octreeArray.resize(levelId[octreeDepth + 1], false);
	cutArray.resize(MAX_NUM);
	cutArray1.resize(MAX_NUM);
	getLevel.resize(levelId[octreeDepth + 1]);
	//for (int i = 0; i < getLevel.size(); i++)
	//	for (int j = octreeDepth; j > -1; j--)
	//		if (levelId[j] <= i) {
	//			getLevel[i] = j; break;
	//		}
	getLevel[0] = 0;
	for (int i = 1; i < octreeDepth; i++)
	{
		std::fill(getLevel.begin() + levelId[i], getLevel.begin() + levelId[i + 1], i);
	}
	ReadRawData(inputFileName, outputFileName);
}

inline void hexGen::ReadRawData(const char* inputFileName, const char* outputFileName) {
	FILE* dataFile = fopen(inputFileName, "r");
	if (NULL == dataFile) {
		std::cerr << "ErrorCode 0: Wrong file name " << inputFileName << std::endl;
		return;
	}
	
	char line[256];
	int i, j, k, l, m, points, elements;
	double box[3][2];// x/y/z min/max
	box[0][0] = MAX_NUM2; box[0][1] = -MAX_NUM2;
	box[1][0] = MAX_NUM2; box[1][1] = -MAX_NUM2;
	box[2][0] = MAX_NUM2; box[2][1] = -MAX_NUM2;
	if (fgets(line, sizeof(line), dataFile) && sscanf(line, "%d %d", &points, &elements) == 2) {
		triMesh.Initialize(elements, points, 3);
		for (i = 0; i < points; i++) triMesh.r[i] = 0;
		triMesh.eNum = elements; triMesh.vNum = points;
		for (i = 0; i < points; i++) {
			fgets(line, sizeof(line), dataFile);
			sscanf(line, "%lf %lf %lf", &triMesh.v[i][0], &triMesh.v[i][1], &triMesh.v[i][2]);
			for (j = 0; j < 3; j++) {
				box[j][0] = (box[j][0] > triMesh.v[i][j]) ? triMesh.v[i][j] : box[j][0];
				box[j][1] = (box[j][1] < triMesh.v[i][j]) ? triMesh.v[i][j] : box[j][1];
			}
		}
		BOX_LENGTH = box[0][1] - box[0][0];
		BOX_LENGTH = (BOX_LENGTH < box[1][1] - box[1][0]) ? box[1][1] - box[1][0] : BOX_LENGTH;
		BOX_LENGTH = (BOX_LENGTH < box[2][1] - box[2][0]) ? box[2][1] - box[2][0] : BOX_LENGTH;
		START_POINT[0] = 0.5 * (box[0][0] + box[0][1] - BOX_LENGTH);
		START_POINT[1] = 0.5 * (box[1][0] + box[1][1] - BOX_LENGTH);
		START_POINT[2] = 0.5 * (box[2][0] + box[2][1] - BOX_LENGTH);
		box[0][0] = MAX_NUM2; box[0][1] = -MAX_NUM2;
		box[1][0] = MAX_NUM2; box[1][1] = -MAX_NUM2;
		box[2][0] = MAX_NUM2; box[2][1] = -MAX_NUM2;
		for (i = 0; i < points; i++) {
			triMesh.v[i][0] -= START_POINT[0];
			triMesh.v[i][1] -= START_POINT[1];
			triMesh.v[i][2] -= START_POINT[2];
			triMesh.v[i][0] *= 100.0 / BOX_LENGTH;
			triMesh.v[i][1] *= 100.0 / BOX_LENGTH;
			triMesh.v[i][2] *= 100.0 / BOX_LENGTH;
			for (j = 0; j < 3; j++) {
				box[j][0] = (box[j][0] > triMesh.v[i][j]) ? triMesh.v[i][j] : box[j][0];
				box[j][1] = (box[j][1] < triMesh.v[i][j]) ? triMesh.v[i][j] : box[j][1];
			}
		}
		for (i = 0; i < elements; i++) {
			fgets(line, sizeof(line), dataFile);
			sscanf(line, "%d %d %d", &triMesh.e[i][0], &triMesh.e[i][1], &triMesh.e[i][2]);
		}
		START_POINT[0] = 0;
		START_POINT[1] = 0;
		START_POINT[2] = 0;
		int publicIdx[2];
		double l1[3], lPublic[3], l2[3], cross1[3], cross2[3], angle;
		for (i = 0; i < triMesh.vNum; i++) triMesh.r[i] = 0;
		for (j = 0; j < triMesh.eNum - 1; j++)
			for (k = 0; k < 3; k++)
				for (l = j + 1; l < triMesh.eNum; l++)
					for (m = 0; m < 3; m++)
						if (triMesh.e[l][m] == triMesh.e[j][k]) {
							if (triMesh.e[j][(k + 1) % 3] == triMesh.e[l][(m + 1) % 3]) {
								publicIdx[0] = 1; publicIdx[1] = 1;
							}
							else if (triMesh.e[j][(k + 1) % 3] == triMesh.e[l][(m + 2) % 3]) {
								publicIdx[0] = 1; publicIdx[1] = 2;
							}
							else if (triMesh.e[j][(k + 2) % 3] == triMesh.e[l][(m + 1) % 3]) {
								publicIdx[0] = 2; publicIdx[1] = 1;
							}
							else if (triMesh.e[j][(k + 2) % 3] == triMesh.e[l][(m + 2) % 3]) {
								publicIdx[0] = 2; publicIdx[1] = 2;
							}
							else continue;
							l1[0] = -triMesh.v[triMesh.e[j][k]][0] + triMesh.v[triMesh.e[j][(k + 3 - publicIdx[0]) % 3]][0];
							l1[1] = -triMesh.v[triMesh.e[j][k]][1] + triMesh.v[triMesh.e[j][(k + 3 - publicIdx[0]) % 3]][1];
							l1[2] = -triMesh.v[triMesh.e[j][k]][2] + triMesh.v[triMesh.e[j][(k + 3 - publicIdx[0]) % 3]][2];
							l2[0] = -triMesh.v[triMesh.e[j][k]][0] + triMesh.v[triMesh.e[l][(m + 3 - publicIdx[1]) % 3]][0];
							l2[1] = -triMesh.v[triMesh.e[j][k]][1] + triMesh.v[triMesh.e[l][(m + 3 - publicIdx[1]) % 3]][1];
							l2[2] = -triMesh.v[triMesh.e[j][k]][2] + triMesh.v[triMesh.e[l][(m + 3 - publicIdx[1]) % 3]][2];
							lPublic[0] = -triMesh.v[triMesh.e[j][k]][0] + triMesh.v[triMesh.e[j][(k + publicIdx[0]) % 3]][0];
							lPublic[1] = -triMesh.v[triMesh.e[j][k]][1] + triMesh.v[triMesh.e[j][(k + publicIdx[0]) % 3]][1];
							lPublic[2] = -triMesh.v[triMesh.e[j][k]][2] + triMesh.v[triMesh.e[j][(k + publicIdx[0]) % 3]][2];
							CROSS(cross1, l1, lPublic)
							CROSS(cross2, l2, lPublic)
							angle = DOT(cross1, cross2) / sqrt(dist(0, 0, 0, cross1) * dist(0, 0, 0, cross2));
							angle = (angle >= -1) ? angle : -1;
							angle = (angle <= 1) ? acos(angle) : 0;
							triMesh.r[triMesh.e[j][k]] += (angle - PI) * (angle - PI);
							break;
						}
		triMesh.WriteToVtk(outputFileName, 1, START_POINT, 5, true);

		BOX_LENGTH = 100;
		START_POINT[0] = 0;
		START_POINT[1] = 0;
		START_POINT[2] = 0;
		BOX_LENGTH_RATIO = BOX_LENGTH / voxelSize;
	}
	else
		std::cerr << "ErrorCode 2: Cannot get total point/element number" << std::endl;
	fclose(dataFile);
}

inline void hexGen::GetCellValue() {
	int i, j, k;
	double center[3], tmp[3], dir[3], len;
	for (i = 0; i < triMesh.eNum; i++) {
		for (j = 0; j < 3; j++)
			// curvature
			if (triMesh.r[triMesh.e[i][j]] > C_THRES[0]) {
				refineTri0.push_back(i);
				refineTriPt0.push_back(triMesh.e[i][j]);
				if (triMesh.r[triMesh.e[i][j]] > C_THRES[1]) {
					refineTri1.push_back(i);
					refineTriPt1.push_back(triMesh.e[i][j]);
					if (triMesh.r[triMesh.e[i][j]] > C_THRES[2]) {
						refineTri2.push_back(i);
						refineTriPt2.push_back(triMesh.e[i][j]);
						if (triMesh.r[triMesh.e[i][j]] > C_THRES[3]) {
							refineTri3.push_back(i);
							refineTriPt3.push_back(triMesh.e[i][j]);
							if (triMesh.r[triMesh.e[i][j]] > C_THRES[4]) {
								refineTri4.push_back(i);
								refineTriPt4.push_back(triMesh.e[i][j]);
								//if (triMesh.r[triMesh.e[i][j]] > C_THRES[5]) {
								//	refineTri5.push_back(i);
								//	refineTriPt5.push_back(triMesh.e[i][j]);
								//}
							}
						}
					}
				}
			}

		// thickness
		center[0] = triMesh.v[triMesh.e[i][1]][0] - triMesh.v[triMesh.e[i][0]][0];
		center[1] = triMesh.v[triMesh.e[i][1]][1] - triMesh.v[triMesh.e[i][0]][1];
		center[2] = triMesh.v[triMesh.e[i][1]][2] - triMesh.v[triMesh.e[i][0]][2];
		tmp[0] = triMesh.v[triMesh.e[i][2]][0] - triMesh.v[triMesh.e[i][0]][0];
		tmp[1] = triMesh.v[triMesh.e[i][2]][1] - triMesh.v[triMesh.e[i][0]][1];
		tmp[2] = triMesh.v[triMesh.e[i][2]][2] - triMesh.v[triMesh.e[i][0]][2];
		CROSS(dir, center, tmp)
			len = sqrt(dist(0, 0, 0, dir));
		dir[0] /= len; dir[1] /= len; dir[2] /= len;
		center[0] = (triMesh.v[triMesh.e[i][0]][0] + triMesh.v[triMesh.e[i][1]][0] + triMesh.v[triMesh.e[i][2]][0]) / 3;
		center[1] = (triMesh.v[triMesh.e[i][0]][1] + triMesh.v[triMesh.e[i][1]][1] + triMesh.v[triMesh.e[i][2]][1]) / 3;
		center[2] = (triMesh.v[triMesh.e[i][0]][2] + triMesh.v[triMesh.e[i][1]][2] + triMesh.v[triMesh.e[i][2]][2]) / 3;

		for (j = i + 1; j < triMesh.eNum; j++) {
			k = Intersect(triMesh.v[triMesh.e[j][0]], triMesh.v[triMesh.e[j][1]], triMesh.v[triMesh.e[j][2]], center, dir, tmp, len);
			tmp[0] = abs(dir[0]) > abs(dir[1]) ? abs(dir[0]) : abs(dir[1]);
			tmp[0] = abs(dir[2]) > tmp[0] ? abs(dir[2]) : tmp[0];
			len = tmp[0] * abs(len);
			if (k == 1)
				if (len < H_THRES[0]) {
					for (k = 0; k < 3; k++) {
						refineTri0.push_back(i); refineTri0.push_back(j); refineTriPt0.push_back(triMesh.e[i][k]);
					}
					if (len < H_THRES[1]) {
						for (k = 0; k < 3; k++) {
							refineTri1.push_back(i); refineTri1.push_back(j); refineTriPt1.push_back(triMesh.e[i][k]);
						}
						if (len < H_THRES[2]) {
							for (k = 0; k < 3; k++) {
								refineTri2.push_back(i); refineTri2.push_back(j); refineTriPt2.push_back(triMesh.e[i][k]);
							}
							if (len < H_THRES[3]) {
								for (k = 0; k < 3; k++) {
									refineTri3.push_back(i); refineTri3.push_back(j); refineTriPt3.push_back(triMesh.e[i][k]);
								}
								if (len < H_THRES[4]) {
									for (k = 0; k < 3; k++) {
										refineTri4.push_back(i); refineTri4.push_back(j); refineTriPt4.push_back(triMesh.e[i][k]);
									}
									//if (len < H_THRES[5]) {
									//	for (k = 0; k < 3; k++) {
									//		refineTri5.push_back(i); refineTri5.push_back(j); refineTriPt5.push_back(triMesh.e[i][k]);
									//	}
									//}
								}
							}
						}
					}
				}
		}
	}
	std::unordered_set<int> s0(refineTriPt0.begin(), refineTriPt0.end());
	refineTriPt0.assign(s0.begin(), s0.end());
	std::unordered_set<int> s1(refineTriPt1.begin(), refineTriPt1.end());
	refineTriPt1.assign(s1.begin(), s1.end());
	std::unordered_set<int> s2(refineTriPt2.begin(), refineTriPt2.end());
	refineTriPt2.assign(s2.begin(), s2.end());
	std::unordered_set<int> s3(refineTriPt3.begin(), refineTriPt3.end());
	refineTriPt3.assign(s3.begin(), s3.end());
	std::unordered_set<int> s4(refineTriPt4.begin(), refineTriPt4.end());
	refineTriPt4.assign(s4.begin(), s4.end());
	//std::unordered_set<int> s5(refineTriPt5.begin(), refineTriPt5.end());
	//refineTriPt5.assign(s5.begin(), s5.end());
	std::unordered_set<int> s00(refineTri0.begin(), refineTri0.end());
	refineTri0.assign(s00.begin(), s00.end());
	std::unordered_set<int> s01(refineTri1.begin(), refineTri1.end());
	refineTri1.assign(s01.begin(), s01.end());
	std::unordered_set<int> s02(refineTri2.begin(), refineTri2.end());
	refineTri2.assign(s02.begin(), s02.end());
	std::unordered_set<int> s03(refineTri3.begin(), refineTri3.end());
	refineTri3.assign(s03.begin(), s03.end());
	std::unordered_set<int> s04(refineTri4.begin(), refineTri4.end());
	refineTri4.assign(s04.begin(), s04.end());
	//std::unordered_set<int> s05(refineTri5.begin(), refineTri5.end());
	//refineTri5.assign(s05.begin(), s05.end());
	for (i = levelId[octreeDepth] - 1; i > -1; i--)
		ComputeCellValue(i, getLevel[i]);
}

inline void hexGen::ComputeCellValue(int octreeId, int level) {
	int i, j, k, l;
	if (level < octreeDepth - 1)
		for (i = 0; i < 8; i++)
			if (octreeArray[Child(octreeId, level, i)]) {
				octreeArray[octreeId] = true;
				for (j = 0; j < 8; j++) {
					k = Child(octreeId, level, j);
					octreeArray[k] = true;
					if (level == octreeDepth - 2 || !octreeArray[Child(k, level + 1, 0)])
						for (l = 0; l < 8; l++) {
							cutArray[leafNum] = Child(k, level + 1, l);
							leafNum++;
						}
				}
				return;
			}

	double box[8][3], center[3], tmp[3], dir[3], len, cellsize = BOX_LENGTH / (1 << level), radius = 0;
	bool flag;

	OctreeidxToXYZ(octreeId, i, j, k, level);

	box[0][0] = START_POINT[0] + (i + 0.5 - CELL_DETECT) * cellsize;
	box[0][1] = START_POINT[1] + (j + 0.5 - CELL_DETECT) * cellsize;
	box[0][2] = START_POINT[2] + (k + 0.5 - CELL_DETECT) * cellsize;
	box[1][0] = START_POINT[0] + (i + 0.5 + CELL_DETECT) * cellsize;
	box[1][1] = START_POINT[1] + (j + 0.5 - CELL_DETECT) * cellsize;
	box[1][2] = START_POINT[2] + (k + 0.5 - CELL_DETECT) * cellsize;
	box[2][0] = START_POINT[0] + (i + 0.5 + CELL_DETECT) * cellsize;
	box[2][1] = START_POINT[1] + (j + 0.5 + CELL_DETECT) * cellsize;
	box[2][2] = START_POINT[2] + (k + 0.5 - CELL_DETECT) * cellsize;
	box[3][0] = START_POINT[0] + (i + 0.5 - CELL_DETECT) * cellsize;
	box[3][1] = START_POINT[1] + (j + 0.5 + CELL_DETECT) * cellsize;
	box[3][2] = START_POINT[2] + (k + 0.5 - CELL_DETECT) * cellsize;
	box[4][0] = START_POINT[0] + (i + 0.5 - CELL_DETECT) * cellsize;
	box[4][1] = START_POINT[1] + (j + 0.5 - CELL_DETECT) * cellsize;
	box[4][2] = START_POINT[2] + (k + 0.5 + CELL_DETECT) * cellsize;
	box[5][0] = START_POINT[0] + (i + 0.5 + CELL_DETECT) * cellsize;
	box[5][1] = START_POINT[1] + (j + 0.5 - CELL_DETECT) * cellsize;
	box[5][2] = START_POINT[2] + (k + 0.5 + CELL_DETECT) * cellsize;
	box[6][0] = START_POINT[0] + (i + 0.5 + CELL_DETECT) * cellsize;
	box[6][1] = START_POINT[1] + (j + 0.5 + CELL_DETECT) * cellsize;
	box[6][2] = START_POINT[2] + (k + 0.5 + CELL_DETECT) * cellsize;
	box[7][0] = START_POINT[0] + (i + 0.5 - CELL_DETECT) * cellsize;
	box[7][1] = START_POINT[1] + (j + 0.5 + CELL_DETECT) * cellsize;
	box[7][2] = START_POINT[2] + (k + 0.5 + CELL_DETECT) * cellsize;

	//if (level == 9) {// 5
	//	for (i = 0; i < refineTriPt5.size(); i++)
	//		if (triMesh.v[refineTriPt5[i]][0] > box[0][0] &&
	//			triMesh.v[refineTriPt5[i]][1] > box[0][1] &&
	//			triMesh.v[refineTriPt5[i]][2] > box[0][2] &&
	//			triMesh.v[refineTriPt5[i]][0] < box[6][0] &&
	//			triMesh.v[refineTriPt5[i]][1] < box[6][1] &&
	//			triMesh.v[refineTriPt5[i]][2] < box[6][2]) {
	//			octreeArray[octreeId] = true; return;
	//		}
	//	for (i = 0; i < refineTri5.size(); i++)
	//		if (Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[0], box[1], box[2]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[0], box[2], box[3]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[4], box[5], box[6]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[4], box[6], box[7]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[0], box[1], box[5]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[0], box[5], box[4]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[1], box[2], box[5]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[2], box[5], box[6]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[2], box[3], box[7]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[2], box[7], box[6]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[0], box[3], box[4]) ||
	//			Intersect(triMesh.v[triMesh.e[refineTri5[i]][0]], triMesh.v[triMesh.e[refineTri5[i]][1]], triMesh.v[triMesh.e[refineTri5[i]][2]], box[3], box[4], box[7])) {
	//			octreeArray[octreeId] = true; return;
	//		}
	//}
	if (level == 8) {// 4
		for (i = 0; i < refineTriPt4.size(); i++)
			if (triMesh.v[refineTriPt4[i]][0] > box[0][0] &&
				triMesh.v[refineTriPt4[i]][1] > box[0][1] &&
				triMesh.v[refineTriPt4[i]][2] > box[0][2] &&
				triMesh.v[refineTriPt4[i]][0] < box[6][0] &&
				triMesh.v[refineTriPt4[i]][1] < box[6][1] &&
				triMesh.v[refineTriPt4[i]][2] < box[6][2]) {
				octreeArray[octreeId] = true; return;
			}
		for (i = 0; i < refineTri4.size(); i++)
			if (Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri4[i]][0]], triMesh.v[triMesh.e[refineTri4[i]][1]], triMesh.v[triMesh.e[refineTri4[i]][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
	}
	else if (level == 7) {// 3
		for (i = 0; i < refineTriPt3.size(); i++)
			if (triMesh.v[refineTriPt3[i]][0] > box[0][0] &&
				triMesh.v[refineTriPt3[i]][1] > box[0][1] &&
				triMesh.v[refineTriPt3[i]][2] > box[0][2] &&
				triMesh.v[refineTriPt3[i]][0] < box[6][0] &&
				triMesh.v[refineTriPt3[i]][1] < box[6][1] &&
				triMesh.v[refineTriPt3[i]][2] < box[6][2]) {
				octreeArray[octreeId] = true; return;
			}
		for (i = 0; i < refineTri3.size(); i++)
			if (Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri3[i]][0]], triMesh.v[triMesh.e[refineTri3[i]][1]], triMesh.v[triMesh.e[refineTri3[i]][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
	}
	else if (level == 6) {// 2
		for (i = 0; i < refineTriPt2.size(); i++)
			if (triMesh.v[refineTriPt2[i]][0] > box[0][0] &&
				triMesh.v[refineTriPt2[i]][1] > box[0][1] &&
				triMesh.v[refineTriPt2[i]][2] > box[0][2] &&
				triMesh.v[refineTriPt2[i]][0] < box[6][0] &&
				triMesh.v[refineTriPt2[i]][1] < box[6][1] &&
				triMesh.v[refineTriPt2[i]][2] < box[6][2]) {
				octreeArray[octreeId] = true; return;
			}
		for (i = 0; i < refineTri2.size(); i++)
			if (Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri2[i]][0]], triMesh.v[triMesh.e[refineTri2[i]][1]], triMesh.v[triMesh.e[refineTri2[i]][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
	}
	else if (level == 5) {// 1
		for (i = 0; i < refineTriPt1.size(); i++)
			if (triMesh.v[refineTriPt1[i]][0] > box[0][0] &&
				triMesh.v[refineTriPt1[i]][1] > box[0][1] &&
				triMesh.v[refineTriPt1[i]][2] > box[0][2] &&
				triMesh.v[refineTriPt1[i]][0] < box[6][0] &&
				triMesh.v[refineTriPt1[i]][1] < box[6][1] &&
				triMesh.v[refineTriPt1[i]][2] < box[6][2]) {
				octreeArray[octreeId] = true; return;
			}
		for (i = 0; i < refineTri1.size(); i++)
			if (Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri1[i]][0]], triMesh.v[triMesh.e[refineTri1[i]][1]], triMesh.v[triMesh.e[refineTri1[i]][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
	}
	else if (level == 4) {// 0
		for (i = 0; i < refineTriPt0.size(); i++)
			if (triMesh.v[refineTriPt0[i]][0] > box[0][0] &&
				triMesh.v[refineTriPt0[i]][1] > box[0][1] &&
				triMesh.v[refineTriPt0[i]][2] > box[0][2] &&
				triMesh.v[refineTriPt0[i]][0] < box[6][0] &&
				triMesh.v[refineTriPt0[i]][1] < box[6][1] &&
				triMesh.v[refineTriPt0[i]][2] < box[6][2]) {
				octreeArray[octreeId] = true; return;
			}
		for (i = 0; i < refineTri0.size(); i++)
			if (Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[refineTri0[i]][0]], triMesh.v[triMesh.e[refineTri0[i]][1]], triMesh.v[triMesh.e[refineTri0[i]][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
	}
	else if (level == 3)// refine all
		for (i = 0; i < triMesh.eNum; i++) {
			for (j = 0; j < 3; j++)
				if (triMesh.v[triMesh.e[i][j]][0] > box[0][0] &&
					triMesh.v[triMesh.e[i][j]][1] > box[0][1] &&
					triMesh.v[triMesh.e[i][j]][2] > box[0][2] &&
					triMesh.v[triMesh.e[i][j]][0] < box[6][0] &&
					triMesh.v[triMesh.e[i][j]][1] < box[6][1] &&
					triMesh.v[triMesh.e[i][j]][2] < box[6][2]) {
					octreeArray[octreeId] = true; return;
				}
			if (Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[0], box[1], box[2]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[0], box[2], box[3]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[4], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[4], box[6], box[7]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[0], box[1], box[5]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[0], box[5], box[4]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[1], box[2], box[5]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[2], box[5], box[6]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[2], box[3], box[7]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[2], box[7], box[6]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[0], box[3], box[4]) ||
				Intersect(triMesh.v[triMesh.e[i][0]], triMesh.v[triMesh.e[i][1]], triMesh.v[triMesh.e[i][2]], box[3], box[4], box[7])) {
				octreeArray[octreeId] = true; return;
			}
		}
}

void hexGen::ConstructOctree() {
	GetCellValue();
	StrongBalancedOctree();
	octreeArray.clear(); cutArray1.clear();
}

inline void hexGen::OctreeidxToXYZ(int octreeId, int& x, int& y, int& z, int level) {
	int idx = octreeId - levelId[level];
	int lRes = levelRes[level];// level of resolution

	x = idx % lRes;
	y = (idx / lRes) % lRes;
	z = idx / (lRes * lRes);
}

inline void hexGen::RefineBrothers(int octreeId, int level, int* octreeIdx) {
	// parent XYZ
	int base = levelId[level];
	int idx = octreeId - base;
	int lRes = levelRes[level], lRes2 = lRes * lRes;
	int x = idx % lRes / 2 * 2;
	int y = (idx / lRes) % lRes / 2 * 2;
	int z = idx / lRes2 / 2 * 2;
	octreeIdx[0] = base + z * lRes2 + y * lRes + x;
	octreeIdx[1] = base + z * lRes2 + y * lRes + x + 1;
	octreeIdx[2] = base + z * lRes2 + (y + 1) * lRes + x;
	octreeIdx[3] = base + z * lRes2 + (y + 1) * lRes + x + 1;
	octreeIdx[4] = base + (z + 1) * lRes2 + y * lRes + x;
	octreeIdx[5] = base + (z + 1) * lRes2 + y * lRes + x + 1;
	octreeIdx[6] = base + (z + 1) * lRes2 + (y + 1) * lRes + x;
	octreeIdx[7] = base + (z + 1) * lRes2 + (y + 1) * lRes + x + 1;

	for (int i = 0; i < 8; i++)
		octreeArray[octreeIdx[i]] = true;
}

inline int hexGen::Child(int octreeId, int level, int i) {
	int idx = octreeId - levelId[level];
	int lRes = levelRes[level], lRes2 = lRes * lRes;
	int x = idx % lRes * 2;
	int y = (idx / lRes) % lRes * 2;
	int z = idx / lRes2 * 2;
	int retIdx;
	int base = levelId[level + 1];
	lRes = levelRes[level + 1]; lRes2 = lRes * lRes;
	switch (i) {
	case 0:
		retIdx = base + z * lRes2 + y * lRes + x;
		break;
	case 1:
		retIdx = base + z * lRes2 + y * lRes + x + 1;
		break;
	case 2:
		retIdx = base + z * lRes2 + (y + 1) * lRes + x;
		break;
	case 3:
		retIdx = base + z * lRes2 + (y + 1) * lRes + x + 1;
		break;
	case 4:
		retIdx = base + (z + 1) * lRes2 + y * lRes + x;
		break;
	case 5:
		retIdx = base + (z + 1) * lRes2 + y * lRes + x + 1;
		break;
	case 6:
		retIdx = base + (z + 1) * lRes2 + (y + 1) * lRes + x;
		break;
	case 7:
		retIdx = base + (z + 1) * lRes2 + (y + 1) * lRes + x + 1;
		break;
	}
	return retIdx;
}

inline void hexGen::StrongBalancedOctree() {
	int unbalancedNode = 0;
	int i, j, k, l;
	int x, y, z;
	int xx, yy, zz;
	int level, level1, level2;
	int leafcellId;
	int cellSize;
	int eightCell[8];
	std::vector<std::vector<std::vector<std::vector<int>>>> preVecEightCell(voxelSize + 1, std::vector<std::vector<std::vector<int>>>(voxelSize + 1, std::vector<std::vector<int>>(voxelSize + 1, std::vector<int>())));

	for (i = 0; i < leafNum; i++) {
		leafcellId = cutArray[i];
		level = getLevel[leafcellId];
		cellSize = voxelSize / (1 << level);
		OctreeidxToXYZ(leafcellId, x, y, z, level);

		for (j = 0; j < 8; j++) {
			switch (j) {
			case 0:
				xx = x * cellSize;
				yy = y * cellSize;
				zz = z * cellSize;
				break;
			case 1:
				xx = (x + 1) * cellSize;
				yy = y * cellSize;
				zz = z * cellSize;
				break;
			case 2:
				xx = (x + 1) * cellSize;
				yy = y * cellSize;
				zz = (z + 1) * cellSize;
				break;
			case 3:
				xx = x * cellSize;
				yy = y * cellSize;
				zz = (z + 1) * cellSize;
				break;
			case 4:
				xx = x * cellSize;
				yy = (y + 1) * cellSize;
				zz = z * cellSize;
				break;
			case 5:
				xx = (x + 1) * cellSize;
				yy = (y + 1) * cellSize;
				zz = z * cellSize;
				break;
			case 6:
				xx = (x + 1) * cellSize;
				yy = (y + 1) * cellSize;
				zz = (z + 1) * cellSize;
				break;
			case 7:
				xx = x * cellSize;
				yy = (y + 1) * cellSize;
				zz = (z + 1) * cellSize;
				break;
			}
			preVecEightCell[xx][yy][zz].push_back(leafcellId);
		}
	}

	for (i = 0; i < leafNum; i++) {
		leafcellId = cutArray[i];
		level = getLevel[leafcellId];
		OctreeidxToXYZ(leafcellId, x, y, z, level);

		cellSize = voxelSize / (1 << level);
		xx = x * cellSize;
		yy = y * cellSize;
		zz = z * cellSize;

		if (preVecEightCell[xx][yy][zz].size() == 8) {
			level1 = 0;
			for (j = 0; j < 8; j++) {
				level2 = getLevel[preVecEightCell[xx][yy][zz][j]];
				level1 = (level1 < level2 ? level2 : level1);
			}

			for (j = 0; j < 8; j++) {
				level2 = getLevel[preVecEightCell[xx][yy][zz][j]];
				if (level1 - level2 > 1 && !octreeArray[preVecEightCell[xx][yy][zz][j]]) {
					octreeArray[preVecEightCell[xx][yy][zz][j]] = true;
					RefineBrothers(preVecEightCell[xx][yy][zz][j], level2, eightCell);
					for (k = 0; k < 8; k++)
						for (l = 0; l < 8; l++) {
							cutArray.push_back(Child(eightCell[k], getLevel[eightCell[k]], l));
							leafNum++;
						}
					unbalancedNode++;
				}
			}
		}
	}

	int leafNum1 = 0;

	for (i = 0; i < cutArray.size(); i++)
		if (!octreeArray[cutArray[i]]) {
			cutArray1[leafNum1] = cutArray[i];
			leafNum1++;
		}

	cutArray = cutArray1;
	leafNum = leafNum1;

	std::cout << "Number of unbalanced nodes: " << unbalancedNode << std::endl;

	if (unbalancedNode > 0)
		StrongBalancedOctree();
}

void hexGen::OutputOctree(const char* fileName) {
	int octreeId;
	int level;
	int cellsize;
	int x, y, z;
	int i, j, k;
	bool overlap;// check if two points overlap with each other
	octreeMesh.Initialize(2 * leafNum);
	octreeENum = 2 * leafNum;
	octreeMesh.eNum = leafNum;

	for (i = 0; i < leafNum; i++)
	{
		octreeId = cutArray[i];
		level = getLevel[octreeId];
		cellsize = voxelSize / (1 << level);
		OctreeidxToXYZ(octreeId, x, y, z, level);
		// cube shape
		//	  7------6
		//	 /|     /|
		//	4-+----5 |
		//	| 3----+-2
		//	|/     |/
		//	0------1

		for (j = 0; j < 8; j++) {
			octreeMesh.v[octreeMesh.vNum][0] = x * cellsize;
			octreeMesh.v[octreeMesh.vNum][1] = y * cellsize;
			octreeMesh.v[octreeMesh.vNum][2] = z * cellsize;
			switch (j) {
			case 1:
				octreeMesh.v[octreeMesh.vNum][0] += cellsize;
				break;
			case 2:
				octreeMesh.v[octreeMesh.vNum][0] += cellsize;
				octreeMesh.v[octreeMesh.vNum][1] += cellsize;
				break;
			case 3:
				octreeMesh.v[octreeMesh.vNum][1] += cellsize;
				break;
			case 4:
				octreeMesh.v[octreeMesh.vNum][2] += cellsize;
				break;
			case 5:
				octreeMesh.v[octreeMesh.vNum][0] += cellsize;
				octreeMesh.v[octreeMesh.vNum][2] += cellsize;
				break;
			case 6:
				octreeMesh.v[octreeMesh.vNum][0] += cellsize;
				octreeMesh.v[octreeMesh.vNum][1] += cellsize;
				octreeMesh.v[octreeMesh.vNum][2] += cellsize;
				break;
			case 7:
				octreeMesh.v[octreeMesh.vNum][1] += cellsize;
				octreeMesh.v[octreeMesh.vNum][2] += cellsize;
				break;
			}
			overlap = false;
			for (k = 0; k < octreeMesh.vNum; k++) {
				if (dist(octreeMesh.v[octreeMesh.vNum], octreeMesh.v[k]) < DIST_THRES) {// overlap
					overlap = true;
					octreeMesh.e[i][j] = k;
					break;
				}
			}
			if (!overlap) {// create new point
				octreeMesh.e[i][j] = octreeMesh.vNum;
				octreeMesh.vNum++;
			}
		}
	}
	cutArray.clear();
	octreeMesh.WriteToVtk(fileName, BOX_LENGTH_RATIO, START_POINT);
}

void hexGen::ReadOctree(const char* inputFileName) {
	FILE* dataFile = fopen(inputFileName, "r");
	if (NULL == dataFile)
	{
		std::cerr << "ErrorCode 0: Wrong file name " << inputFileName << std::endl;
		return;
	}

	char line[256];
	int i, j;
	for (i = 0; i < 4 && fgets(line, sizeof(line), dataFile); i++) {
	}

	int points;
	if (fgets(line, sizeof(line), dataFile) && sscanf(line, "POINTS %d double", &points) == 1) {
		octreeMesh.Initialize(points);
		octreeENum = points;
		octreeMesh.vNum = points;
		for (i = 0; i < points; i++) {
			fgets(line, sizeof(line), dataFile);
			sscanf(line, "%lf %lf %lf", &octreeMesh.v[i][0], &octreeMesh.v[i][1], &octreeMesh.v[i][2]);
			octreeMesh.v[i][0] = round((octreeMesh.v[i][0] - START_POINT[0]) / BOX_LENGTH_RATIO);
			octreeMesh.v[i][1] = round((octreeMesh.v[i][1] - START_POINT[1]) / BOX_LENGTH_RATIO);
			octreeMesh.v[i][2] = round((octreeMesh.v[i][2] - START_POINT[2]) / BOX_LENGTH_RATIO);
		}
		if (fgets(line, sizeof(line), dataFile) && sscanf(line, "CELLS %d %d", &leafNum, &i) == 2) {
			octreeMesh.eNum = leafNum;
			for (i = 0; i < leafNum; i++) {
				fgets(line, sizeof(line), dataFile);
				sscanf(line, "%d %d %d %d %d %d %d %d %d", &j, &octreeMesh.e[i][0], &octreeMesh.e[i][1], &octreeMesh.e[i][2], &octreeMesh.e[i][3], &octreeMesh.e[i][4], &octreeMesh.e[i][5], &octreeMesh.e[i][6], &octreeMesh.e[i][7]);
			}
		}
		else
			std::cerr << "ErrorCode 2: Cannot get total element number" << std::endl;
	}
	else
		std::cerr << "ErrorCode 1: Cannot get total point number" << std::endl;
	fclose(dataFile);
}

void hexGen::DualFullHexMeshExtraction(const char* fileName) {
	InitiateElementValence();

	hexMesh.Initialize(2 * leafNum);
	hexMeshENum = 2 * leafNum;
	int i, j, k, l = 0, m, vNum;
	double posI, posJ, delta, p[32][3], p16[16][3], ptmp[3], z[3];
	bool overlap, stepI, stepJ, checkOverlap;

	// tmp: center of all octree cells
	double(*tmp)[3];
	tmp = new double[leafNum][3];
	for (i = 0; i < leafNum; i++) {
		tmp[i][0] = 0.5 * (octreeMesh.v[octreeMesh.e[i][0]][0] + octreeMesh.v[octreeMesh.e[i][6]][0]);
		tmp[i][1] = 0.5 * (octreeMesh.v[octreeMesh.e[i][0]][1] + octreeMesh.v[octreeMesh.e[i][6]][1]);
		tmp[i][2] = 0.5 * (octreeMesh.v[octreeMesh.e[i][0]][2] + octreeMesh.v[octreeMesh.e[i][6]][2]);
	}

	// empty hex candidates + generate regular hex
	int countValence, num[8], collectNumLength = 0; double size[8];
	std::vector<std::vector<int>> collectNum(octreeMesh.vNum, std::vector<int>(9));
	for (i = 0; i < octreeMesh.vNum; i++) {
		countValence = 0;
		for (j = 0; j < leafNum; j++)
			for (k = 0; k < 8; k++)
				if (octreeMesh.e[j][k] == i) {
					size[countValence] = abs(octreeMesh.v[octreeMesh.e[j][0]][0] - octreeMesh.v[octreeMesh.e[j][6]][0]);
					num[idTransform[k]] = j;
					countValence++;
				}
		if (countValence == 8 && !(size[0] == size[1] && size[0] == size[2] && size[0] == size[3] && size[0] == size[4] && size[0] == size[5] && size[0] == size[6] && size[0] == size[7])) {
			collectNum[collectNumLength][8] = i;
			for (m = 0; m < 8; m++)
				collectNum[collectNumLength][m] = num[m];
			collectNumLength += 1;
		}
	}

	// convert all templates to hexahedral elements
	for (i = 0; i < leafNum; i++) {
		if (elementValenceNumber[i][3] == 1 && elementValenceNumber[elementValence[i][3][0]][2] == 1 &&
			elementValenceNumber[i][4] == 1 && elementValenceNumber[elementValence[i][4][0]][1] == 1 &&
			elementValenceNumber[i][5] == 1 && elementValenceNumber[elementValence[i][5][0]][0] == 1 &&
			elementValenceNumber[elementValence[i][3][0]][5] == 1 && elementValenceNumber[elementValence[elementValence[i][3][0]][5][0]][0] == 1 &&
			elementValenceNumber[elementValence[i][4][0]][5] == 1 && elementValenceNumber[elementValence[elementValence[i][4][0]][5][0]][0] == 1 &&
			elementValenceNumber[elementValence[i][3][0]][4] == 1 && elementValenceNumber[elementValence[elementValence[i][3][0]][4][0]][1] == 1 &&
			elementValenceNumber[elementValence[elementValence[i][3][0]][4][0]][5] == 1 && elementValenceNumber[elementValence[elementValence[elementValence[i][3][0]][4][0]][5][0]][0] == 1) {// construct 1 element
			p[0][0] = tmp[i][0]; p[0][1] = tmp[i][1]; p[0][2] = tmp[i][2];
			p[1][0] = tmp[elementValence[i][3][0]][0]; p[1][1] = tmp[elementValence[i][3][0]][1]; p[1][2] = tmp[elementValence[i][3][0]][2];
			p[2][0] = tmp[elementValence[elementValence[i][3][0]][4][0]][0]; p[2][1] = tmp[elementValence[elementValence[i][3][0]][4][0]][1]; p[2][2] = tmp[elementValence[elementValence[i][3][0]][4][0]][2];
			p[3][0] = tmp[elementValence[i][4][0]][0]; p[3][1] = tmp[elementValence[i][4][0]][1]; p[3][2] = tmp[elementValence[i][4][0]][2];
			p[4][0] = tmp[elementValence[i][5][0]][0]; p[4][1] = tmp[elementValence[i][5][0]][1]; p[4][2] = tmp[elementValence[i][5][0]][2];
			p[5][0] = tmp[elementValence[elementValence[i][3][0]][5][0]][0]; p[5][1] = tmp[elementValence[elementValence[i][3][0]][5][0]][1]; p[5][2] = tmp[elementValence[elementValence[i][3][0]][5][0]][2];
			p[6][0] = tmp[elementValence[elementValence[elementValence[i][3][0]][5][0]][4][0]][0]; p[6][1] = tmp[elementValence[elementValence[elementValence[i][3][0]][5][0]][4][0]][1]; p[6][2] = tmp[elementValence[elementValence[elementValence[i][3][0]][5][0]][4][0]][2];
			p[7][0] = tmp[elementValence[elementValence[i][4][0]][5][0]][0]; p[7][1] = tmp[elementValence[elementValence[i][4][0]][5][0]][1]; p[7][2] = tmp[elementValence[elementValence[i][4][0]][5][0]][2];
			for (j = 0; j < 8; j++) {
				overlap = false;
				for (k = 0; k < hexMesh.vNum; k++)
					if (dist(hexMesh.v[k], p[j]) < DIST_THRES) {// overlap
						overlap = true;
						hexMesh.e[hexMesh.eNum][j] = k;
						break;
					}
				if (!overlap) {// create new point
					hexMesh.e[hexMesh.eNum][j] = hexMesh.vNum;
					hexMesh.v[hexMesh.vNum][0] = p[j][0];
					hexMesh.v[hexMesh.vNum][1] = p[j][1];
					hexMesh.v[hexMesh.vNum][2] = p[j][2];
					hexMesh.vNum++;
				}
			}
			hexMesh.eNum++;
		}
		for (j = 0; j < 6; j++)
			if (elementValenceNumber[i][j] == 4) {// construct 13 element template
				// the template can be at four possible positions, need to check where it is with pos
				// the checking is very tricky
				// reduce tmp[i][0-3] in both directions with the size of the big cube, calculate the step it takes to get to below 0: stepI and stepJ
				// if stepI is true, the position will be in 0/3 otherwise in 1/2
				// if stepJ is true, the position will be in 0/1 otherwise in 2/3
				stepI = false; stepJ = false; delta = octreeMesh.v[octreeMesh.e[i][1]][0] - octreeMesh.v[octreeMesh.e[i][0]][0];
				if (j == 0 || j == 5) {// z direction
					posI = tmp[i][0]; posJ = tmp[i][1];
				}
				else if (j == 1 || j == 4) {// y direction
					posI = tmp[i][0]; posJ = tmp[i][2];
				}
				else {// x direction
					posI = tmp[i][1]; posJ = tmp[i][2];
				}
				while (posI > 0) {
					posI -= delta; stepI = !stepI;
				}
				while (posJ > 0) {
					posJ -= delta; stepJ = !stepJ;
				}
				if (stepI && stepJ) {
					// there are 32 vertices in total
					// 12      13/30    14/31    15
					//      28                29
					// 8/24    9/25     10/26    11/27
					//
					// 4/20	   5/21     6/22     7/23
					//      18                19
					// 0       1/16     2/17     3
					// there are 13 hexahedrals in total
					// 0:	0	1	5	4	18	16	21	20
					// 1:	1	2	6	5	16	17	22	21
					// 2:	2	3	7	6	17	19	23	22
					// 3:	4	5	9	8	20	21	25	24
					// 4:	5	6	10	9	21	22	26	25
					// 5:	6	7	11	10	22	23	27	26
					// 6:	8	9	13	12	24	25	30	28
					// 7:	9	10	14	13	25	26	31	30
					// 8:	10	11	15	14	26	27	29	31
					// 9:	20	21	25	24	18	16	30	28
					// 10:	22	23	27	26	17	19	29	31
					// 11:	21	22	26	25	16	17	31	30
					// 12:	16	17	31	30	18	19	29	28
					// points
					p[0][0] = tmp[elementValence[i][j][0]][0]; p[0][1] = tmp[elementValence[i][j][0]][1]; p[0][2] = tmp[elementValence[i][j][0]][2];
					p[1][0] = tmp[elementValence[i][j][1]][0]; p[1][1] = tmp[elementValence[i][j][1]][1]; p[1][2] = tmp[elementValence[i][j][1]][2];
					p[2][0] = 2 * p[1][0] - p[0][0]; p[2][1] = 2 * p[1][1] - p[0][1]; p[2][2] = 2 * p[1][2] - p[0][2];
					p[3][0] = 2 * p[2][0] - p[1][0]; p[3][1] = 2 * p[2][1] - p[1][1]; p[3][2] = 2 * p[2][2] - p[1][2];
					p[4][0] = tmp[elementValence[i][j][3]][0]; p[4][1] = tmp[elementValence[i][j][3]][1]; p[4][2] = tmp[elementValence[i][j][3]][2];
					p[5][0] = tmp[elementValence[i][j][2]][0]; p[5][1] = tmp[elementValence[i][j][2]][1]; p[5][2] = tmp[elementValence[i][j][2]][2];
					p[6][0] = 2 * p[5][0] - p[4][0]; p[6][1] = 2 * p[5][1] - p[4][1]; p[6][2] = 2 * p[5][2] - p[4][2];
					p[7][0] = 2 * p[6][0] - p[5][0]; p[7][1] = 2 * p[6][1] - p[5][1]; p[7][2] = 2 * p[6][2] - p[5][2];
					p[8][0] = 2 * p[4][0] - p[0][0]; p[8][1] = 2 * p[4][1] - p[0][1]; p[8][2] = 2 * p[4][2] - p[0][2];
					p[9][0] = 2 * p[5][0] - p[1][0]; p[9][1] = 2 * p[5][1] - p[1][1]; p[9][2] = 2 * p[5][2] - p[1][2];
					p[10][0] = 2 * p[5][0] - p[0][0]; p[10][1] = 2 * p[5][1] - p[0][1]; p[10][2] = 2 * p[5][2] - p[0][2];
					p[11][0] = 2 * p[10][0] - p[9][0]; p[11][1] = 2 * p[10][1] - p[9][1]; p[11][2] = 2 * p[10][2] - p[9][2];
					p[12][0] = 2 * p[8][0] - p[4][0]; p[12][1] = 2 * p[8][1] - p[4][1]; p[12][2] = 2 * p[8][2] - p[4][2];
					p[13][0] = 2 * p[9][0] - p[5][0]; p[13][1] = 2 * p[9][1] - p[5][1]; p[13][2] = 2 * p[9][2] - p[5][2];
					p[14][0] = 2 * p[13][0] - p[12][0]; p[14][1] = 2 * p[13][1] - p[12][1]; p[14][2] = 2 * p[13][2] - p[12][2];
					p[15][0] = 2 * p[14][0] - p[13][0]; p[15][1] = 2 * p[14][1] - p[13][1]; p[15][2] = 2 * p[14][2] - p[13][2];
					// z points upward
					z[0] = 0.25 * (p[0][0] + p[1][0] + p[4][0] + p[5][0]); z[0] = tmp[i][0] - z[0];
					z[1] = 0.25 * (p[0][1] + p[1][1] + p[4][1] + p[5][1]); z[1] = tmp[i][1] - z[1];
					z[2] = 0.25 * (p[0][2] + p[1][2] + p[4][2] + p[5][2]); z[2] = tmp[i][2] - z[2];
					p[16][0] = p[1][0] + 268 * z[0] / 375 + (p[4][0] - p[0][0]) * 0.072;
					p[16][1] = p[1][1] + 268 * z[1] / 375 + (p[4][1] - p[0][1]) * 0.072;
					p[16][2] = p[1][2] + 268 * z[2] / 375 + (p[4][2] - p[0][2]) * 0.072;
					p[17][0] = p[2][0] + 268 * z[0] / 375 + (p[4][0] - p[0][0]) * 0.072;
					p[17][1] = p[2][1] + 268 * z[1] / 375 + (p[4][1] - p[0][1]) * 0.072;
					p[17][2] = p[2][2] + 268 * z[2] / 375 + (p[4][2] - p[0][2]) * 0.072;
					p[18][0] = tmp[i][0]; p[18][1] = tmp[i][1]; p[18][2] = tmp[i][2];
					p[19][0] = p[18][0] + 2 * (p[1][0] - p[0][0]); p[19][1] = p[18][1] + 2 * (p[1][1] - p[0][1]); p[19][2] = p[18][2] + 2 * (p[1][2] - p[0][2]);
					p[20][0] = p[4][0] + 268 * z[0] / 375 + (p[1][0] - p[0][0]) * 0.072;
					p[20][1] = p[4][1] + 268 * z[1] / 375 + (p[1][1] - p[0][1]) * 0.072;
					p[20][2] = p[4][2] + 268 * z[2] / 375 + (p[1][2] - p[0][2]) * 0.072;
					p[21][0] = p[5][0] + z[0] / 5 - (p[1][0] - p[0][0]) * 0.112 + (p[4][0] - p[0][0]) * 0.056;
					p[21][1] = p[5][1] + z[1] / 5 - (p[1][1] - p[0][1]) * 0.112 + (p[4][1] - p[0][1]) * 0.056;
					p[21][2] = p[5][2] + z[2] / 5 - (p[1][2] - p[0][2]) * 0.112 + (p[4][2] - p[0][2]) * 0.056;
					p[22][0] = p[6][0] + z[0] / 5 + (p[1][0] - p[0][0]) * 0.112 + (p[4][0] - p[0][0]) * 0.056;
					p[22][1] = p[6][1] + z[1] / 5 + (p[1][1] - p[0][1]) * 0.112 + (p[4][1] - p[0][1]) * 0.056;
					p[22][2] = p[6][2] + z[2] / 5 + (p[1][2] - p[0][2]) * 0.112 + (p[4][2] - p[0][2]) * 0.056;
					p[23][0] = p[7][0] + 268 * z[0] / 375 - (p[1][0] - p[0][0]) * 0.072;
					p[23][1] = p[7][1] + 268 * z[1] / 375 - (p[1][1] - p[0][1]) * 0.072;
					p[23][2] = p[7][2] + 268 * z[2] / 375 - (p[1][2] - p[0][2]) * 0.072;
					p[24][0] = p[8][0] + 268 * z[0] / 375 + (p[1][0] - p[0][0]) * 0.072;
					p[24][1] = p[8][1] + 268 * z[1] / 375 + (p[1][1] - p[0][1]) * 0.072;
					p[24][2] = p[8][2] + 268 * z[2] / 375 + (p[1][2] - p[0][2]) * 0.072;
					p[25][0] = p[9][0] + z[0] / 5 - (p[1][0] - p[0][0]) * 0.112 - (p[4][0] - p[0][0]) * 0.056;
					p[25][1] = p[9][1] + z[1] / 5 - (p[1][1] - p[0][1]) * 0.112 - (p[4][1] - p[0][1]) * 0.056;
					p[25][2] = p[9][2] + z[2] / 5 - (p[1][2] - p[0][2]) * 0.112 - (p[4][2] - p[0][2]) * 0.056;
					p[26][0] = p[10][0] + z[0] / 5 + (p[1][0] - p[0][0]) * 0.112 - (p[4][0] - p[0][0]) * 0.056;
					p[26][1] = p[10][1] + z[1] / 5 + (p[1][1] - p[0][1]) * 0.112 - (p[4][1] - p[0][1]) * 0.056;
					p[26][2] = p[10][2] + z[2] / 5 + (p[1][2] - p[0][2]) * 0.112 - (p[4][2] - p[0][2]) * 0.056;
					p[27][0] = p[11][0] + 268 * z[0] / 375 - (p[1][0] - p[0][0]) * 0.072;
					p[27][1] = p[11][1] + 268 * z[1] / 375 - (p[1][1] - p[0][1]) * 0.072;
					p[27][2] = p[11][2] + 268 * z[2] / 375 - (p[1][2] - p[0][2]) * 0.072;
					p[28][0] = p[18][0] + 2 * (p[4][0] - p[0][0]); p[28][1] = p[18][1] + 2 * (p[4][1] - p[0][1]); p[28][2] = p[18][2] + 2 * (p[4][2] - p[0][2]);
					p[29][0] = p[18][0] + 2 * (p[5][0] - p[0][0]); p[29][1] = p[18][1] + 2 * (p[5][1] - p[0][1]); p[29][2] = p[18][2] + 2 * (p[5][2] - p[0][2]);
					p[30][0] = p[13][0] + 268 * z[0] / 375 - (p[4][0] - p[0][0]) * 0.072;
					p[30][1] = p[13][1] + 268 * z[1] / 375 - (p[4][1] - p[0][1]) * 0.072;
					p[30][2] = p[13][2] + 268 * z[2] / 375 - (p[4][2] - p[0][2]) * 0.072;
					p[31][0] = p[14][0] + 268 * z[0] / 375 - (p[4][0] - p[0][0]) * 0.072;
					p[31][1] = p[14][1] + 268 * z[1] / 375 - (p[4][1] - p[0][1]) * 0.072;
					p[31][2] = p[14][2] + 268 * z[2] / 375 - (p[4][2] - p[0][2]) * 0.072;
					// elements
					stepI = j == 1 || j == 3 || j == 5;
					for (k = 0; k < 13; k++) {
						for (l = 0; l < 8; l++) {
							overlap = false;
							for (m = 0; m < hexMesh.vNum; m++)
								if (dist(hexMesh.v[m], p[t1Id[stepI][k][l]]) < DIST_THRES) {// overlap
									overlap = true;
									hexMesh.e[hexMesh.eNum][l] = m;
									break;
								}
							if (!overlap) {// create new point
								hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
								hexMesh.v[hexMesh.vNum][0] = p[t1Id[stepI][k][l]][0];
								hexMesh.v[hexMesh.vNum][1] = p[t1Id[stepI][k][l]][1];
								hexMesh.v[hexMesh.vNum][2] = p[t1Id[stepI][k][l]][2];
								hexMesh.vNum++;
							}
						}
						hexMesh.eNum++;
					}

					// delete corresponding point in collectNum[][8]
					ptmp[0] = 0.5 * (p[21][0] + p[26][0] + z[0] * 4 / 15); ptmp[1] = 0.5 * (p[21][1] + p[26][1] + z[1] * 4 / 15); ptmp[2] = 0.5 * (p[21][2] + p[26][2] + z[2] * 4 / 15);
					for (k = 0; k < collectNumLength; k++)
						if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
							collectNum.erase(collectNum.begin() + k);
							collectNumLength--;
							break;
						}

					// up to now, the 13 element template is successfully constructed
					// need to fill the 4, 3, 5 element transition templates

					if (elementValenceNumber[i][pSId[j][0]] == 1 && elementValenceNumber[elementValence[i][pSId[j][0]][0]][j] == 4) {// 4 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						//    7n    12
						// 6n           28
						//    4n/5n 8/24

						//    2n/3n 4/20
						// 1n           18
						//    0n    0
						// there are 4 hexahedrals in total
						// 0:	0n	0	4	2n	1n	18	20	3n
						// 1:	2n	4	8	4n	3n	20	24	5n
						// 2:	3n	20	24	5n	1n	18	28	6n
						// 3:	4n	8	12	7n	5n	24	28	6n
						// old points: 0, 4, 8, 12, 18, 20, 24, 28 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	8n	9n	2n	1n	12n	13n	3n
						// 1:	2n	9n	10n	4n	3n	13n	14n	5n
						// 2:	3n	13n	14n	5n	1n	12n	15n	6n
						// 3:	4n	10n	11n	7n	5n	14n	15n	6n
						// points
						p16[0][0] = 2 * p[0][0] - p[1][0]; p16[0][1] = 2 * p[0][1] - p[1][1]; p16[0][2] = 2 * p[0][2] - p[1][2];
						p16[1][0] = 2 * p[18][0] - p[19][0]; p16[1][1] = 2 * p[18][1] - p[19][1]; p16[1][2] = 2 * p[18][2] - p[19][2];
						p16[2][0] = 2 * p[4][0] - p[5][0]; p16[2][1] = 2 * p[4][1] - p[5][1]; p16[2][2] = 2 * p[4][2] - p[5][2];
						p16[3][0] = p[20][0] + 1.144 * (p[0][0] - p[1][0]); p16[3][1] = p[20][1] + 1.144 * (p[0][1] - p[1][1]); p16[3][2] = p[20][2] + 1.144 * (p[0][2] - p[1][2]);
						p16[4][0] = 2 * p[8][0] - p[9][0]; p16[4][1] = 2 * p[8][1] - p[9][1]; p16[4][2] = 2 * p[8][2] - p[9][2];
						p16[5][0] = p[24][0] + 1.144 * (p[0][0] - p[1][0]); p16[5][1] = p[24][1] + 1.144 * (p[0][1] - p[1][1]); p16[5][2] = p[24][2] + 1.144 * (p[0][2] - p[1][2]);
						p16[6][0] = 2 * p[28][0] - p[29][0]; p16[6][1] = 2 * p[28][1] - p[29][1]; p16[6][2] = 2 * p[28][2] - p[29][2];
						p16[7][0] = 2 * p[12][0] - p[13][0]; p16[7][1] = 2 * p[12][1] - p[13][1]; p16[7][2] = 2 * p[12][2] - p[13][2];
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[4][0]; p16[9][1] = p[4][1]; p16[9][2] = p[4][2];
						p16[10][0] = p[8][0]; p16[10][1] = p[8][1]; p16[10][2] = p[8][2];
						p16[11][0] = p[12][0]; p16[11][1] = p[12][1]; p16[11][2] = p[12][2];
						p16[12][0] = p[18][0]; p16[12][1] = p[18][1]; p16[12][2] = p[18][2];
						p16[13][0] = p[20][0]; p16[13][1] = p[20][1]; p16[13][2] = p[20][2];
						p16[14][0] = p[24][0]; p16[14][1] = p[24][1]; p16[14][2] = p[24][2];
						p16[15][0] = p[28][0]; p16[15][1] = p[28][1]; p16[15][2] = p[28][2];
						// elements
						for (k = 0; k < 4; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t2Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t2Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t2Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t2Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						ptmp[0] = 0.5 * (p16[2][0] + p[8][0]) + z[0] / 3; ptmp[1] = 0.5 * (p16[2][1] + p[8][1]) + z[1] / 3; ptmp[2] = 0.5 * (p16[2][2] + p[8][2]) + z[2] / 3;
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
					}

					if (elementValenceNumber[i][pSId[j][1]] == 1 && elementValenceNumber[elementValence[i][pSId[j][1]][0]][j] == 4) {// 4 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						//   18             19
						// 0    1/16  2/17     3
						// 0n   2n/3n 4n/5n    7n
						//   1n             6n
						// there are 4 hexahedrals in total
						// 0:	0	0n	2n	1	18	1n	3n	16
						// 1:	1	2n	4n	2	16	3n	5n	17
						// 2:	16	3n	5n	17	18	1n	6n	19
						// 3:	2	4n	7n	3	17	5n	6n	19
						// old points: 0, 1, 2, 3, 16, 17, 18, 19 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	8n	0n	2n	9n	14n	1n	3n	12n
						// 1:	9n	2n	4n	10n	12n	3n	5n	13n
						// 2:	12n	3n	5n	13n	14n	1n	6n	15n
						// 3:	10n	4n	7n	11n	13n	5n	6n	15n
						// points
						p16[0][0] = 2 * p[0][0] - p[4][0]; p16[0][1] = 2 * p[0][1] - p[4][1]; p16[0][2] = 2 * p[0][2] - p[4][2];
						p16[1][0] = 2 * p[18][0] - p[28][0]; p16[1][1] = 2 * p[18][1] - p[28][1]; p16[1][2] = 2 * p[18][2] - p[28][2];
						p16[2][0] = 2 * p[1][0] - p[5][0]; p16[2][1] = 2 * p[1][1] - p[5][1]; p16[2][2] = 2 * p[1][2] - p[5][2];
						p16[3][0] = p[16][0] + 1.144 * (p[0][0] - p[4][0]); p16[3][1] = p[16][1] + 1.144 * (p[0][1] - p[4][1]); p16[3][2] = p[16][2] + 1.144 * (p[0][2] - p[4][2]);
						p16[4][0] = 2 * p[2][0] - p[6][0]; p16[4][1] = 2 * p[2][1] - p[6][1]; p16[4][2] = 2 * p[2][2] - p[6][2];
						p16[5][0] = p[17][0] + 1.144 * (p[0][0] - p[4][0]); p16[5][1] = p[17][1] + 1.144 * (p[0][1] - p[4][1]); p16[5][2] = p[17][2] + 1.144 * (p[0][2] - p[4][2]);
						p16[6][0] = 2 * p[19][0] - p[29][0]; p16[6][1] = 2 * p[19][1] - p[29][1]; p16[6][2] = 2 * p[19][2] - p[29][2];
						p16[7][0] = 2 * p[3][0] - p[7][0]; p16[7][1] = 2 * p[3][1] - p[7][1]; p16[7][2] = 2 * p[3][2] - p[7][2];
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[1][0]; p16[9][1] = p[1][1]; p16[9][2] = p[1][2];
						p16[10][0] = p[2][0]; p16[10][1] = p[2][1]; p16[10][2] = p[2][2];
						p16[11][0] = p[3][0]; p16[11][1] = p[3][1]; p16[11][2] = p[3][2];
						p16[12][0] = p[16][0]; p16[12][1] = p[16][1]; p16[12][2] = p[16][2];
						p16[13][0] = p[17][0]; p16[13][1] = p[17][1]; p16[13][2] = p[17][2];
						p16[14][0] = p[18][0]; p16[14][1] = p[18][1]; p16[14][2] = p[18][2];
						p16[15][0] = p[19][0]; p16[15][1] = p[19][1]; p16[15][2] = p[19][2];
						// elements
						for (k = 0; k < 4; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t22Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t22Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t22Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t22Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						ptmp[0] = 0.5 * (p16[4][0] + p[1][0]) + z[0] / 3; ptmp[1] = 0.5 * (p16[4][1] + p[1][1]) + z[1] / 3; ptmp[2] = 0.5 * (p16[4][2] + p[1][2]) + z[2] / 3;
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
					}

					if (elementValenceNumber[i][pSId[j][0]] == 4) {// 3 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 6n/7n 12/28
						// 4n/5n 8/24
						// 2n/3n 4/20
						// 0n/1n 0/18
						// there are 3 hexahedrals in total
						// 0:	0n	0	4	2n	1n	18	20	3n
						// 1:	2n	4	8	4n	3n	20	24	5n
						// 2:	4n	8	12	6n	5n	24	28	7n
						// old points: 0, 4, 8, 12, 18, 20, 24, 28 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	8n	9n	2n	1n	12n	13n	3n
						// 1:	2n	9n	10n	4n	3n	13n	14n	5n
						// 2:	4n	10n	11n	6n	5n	14n	15n	7n
						// points
						p16[0][0] = 2 * p[0][0] - p[1][0]; p16[0][1] = 2 * p[0][1] - p[1][1]; p16[0][2] = 2 * p[0][2] - p[1][2];
						p16[1][0] = p16[0][0] + 2 * z[0] / 3; p16[1][1] = p16[0][1] + 2 * z[1] / 3; p16[1][2] = p16[0][2] + 2 * z[2] / 3;
						p16[2][0] = 2 * p[4][0] - p[5][0]; p16[2][1] = 2 * p[4][1] - p[5][1]; p16[2][2] = 2 * p[4][2] - p[5][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[8][0] - p[9][0]; p16[4][1] = 2 * p[8][1] - p[9][1]; p16[4][2] = 2 * p[8][2] - p[9][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[6][0] = 2 * p[12][0] - p[13][0]; p16[6][1] = 2 * p[12][1] - p[13][1]; p16[6][2] = 2 * p[12][2] - p[13][2];
						p16[7][0] = p16[6][0] + 2 * z[0] / 3; p16[7][1] = p16[6][1] + 2 * z[1] / 3; p16[7][2] = p16[6][2] + 2 * z[2] / 3;
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[4][0]; p16[9][1] = p[4][1]; p16[9][2] = p[4][2];
						p16[10][0] = p[8][0]; p16[10][1] = p[8][1]; p16[10][2] = p[8][2];
						p16[11][0] = p[12][0]; p16[11][1] = p[12][1]; p16[11][2] = p[12][2];
						p16[12][0] = p[18][0]; p16[12][1] = p[18][1]; p16[12][2] = p[18][2];
						p16[13][0] = p[20][0]; p16[13][1] = p[20][1]; p16[13][2] = p[20][2];
						p16[14][0] = p[24][0]; p16[14][1] = p[24][1]; p16[14][2] = p[24][2];
						p16[15][0] = p[28][0]; p16[15][1] = p[28][1]; p16[15][2] = p[28][2];
						if (elementValenceNumber[elementValence[elementValence[i][pSId[j][0]][0]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][0]][1]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][0]][2]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][0]][3]][j][0]][5 - j] == 4) {
							p16[0][0] = 2 * p[18][0] - p[19][0] - 4 * z[0] / 3; p16[0][1] = 2 * p[18][1] - p[19][1] - 4 * z[1] / 3; p16[0][2] = 2 * p[18][2] - p[19][2] - 4 * z[2] / 3;
							p16[6][0] = 2 * p[28][0] - p[29][0] - 4 * z[0] / 3; p16[6][1] = 2 * p[28][1] - p[29][1] - 4 * z[1] / 3; p16[6][2] = 2 * p[28][2] - p[29][2] - 4 * z[2] / 3;
							p16[2][0] = p16[5][0] + p[4][0] - p[24][0]; p16[2][1] = p16[5][1] + p[4][1] - p[24][1]; p16[2][2] = p16[5][2] + p[4][2] - p[24][2];
							p16[4][0] = p16[5][0] + p[4][0] - p[20][0]; p16[4][1] = p16[5][1] + p[4][1] - p[20][1]; p16[4][2] = p16[5][2] + p[4][2] - p[20][2];
						}
						// elements
						for (k = 0; k < 3; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t3Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t3Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t3Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t3Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						checkOverlap = true;
						ptmp[0] = 0.5 * (p16[5][0] + p[4][0]); ptmp[1] = 0.5 * (p16[5][1] + p[4][1]); ptmp[2] = 0.5 * (p16[5][2] + p[4][2]);
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								checkOverlap = false;
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
						if (checkOverlap)// the template is already created
							hexMesh.eNum -= 3;
					}

					if (elementValenceNumber[i][pSId[j][1]] == 4) {// 3 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 0/18  1/16  2/17  3/19
						// 0n/1n 2n/3n 4n/5n 6n/7n
						// there are 3 hexahedrals in total
						// 0:	0n	2n	1	0	1n	3n	16	18
						// 1:	2n	4n	2	1	3n	5n	17	16
						// 2:	4n	6n	3	2	5n	7n	19	17
						// old points: 0, 1, 2, 3, 16, 17, 18, 19 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	2n	9n	8n	1n	3n	12n	14n
						// 1:	2n	4n	10n	9n	3n	5n	13n	12n
						// 2:	4n	6n	11n	10n	5n	7n	15n	13n
						// points
						p16[0][0] = 2 * p[0][0] - p[4][0]; p16[0][1] = 2 * p[0][1] - p[4][1]; p16[0][2] = 2 * p[0][2] - p[4][2];
						p16[1][0] = p16[0][0] + 2 * z[0] / 3; p16[1][1] = p16[0][1] + 2 * z[1] / 3; p16[1][2] = p16[0][2] + 2 * z[2] / 3;
						p16[2][0] = 2 * p[1][0] - p[5][0]; p16[2][1] = 2 * p[1][1] - p[5][1]; p16[2][2] = 2 * p[1][2] - p[5][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[2][0] - p[6][0]; p16[4][1] = 2 * p[2][1] - p[6][1]; p16[4][2] = 2 * p[2][2] - p[6][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[6][0] = 2 * p[3][0] - p[7][0]; p16[6][1] = 2 * p[3][1] - p[7][1]; p16[6][2] = 2 * p[3][2] - p[7][2];
						p16[7][0] = p16[6][0] + 2 * z[0] / 3; p16[7][1] = p16[6][1] + 2 * z[1] / 3; p16[7][2] = p16[6][2] + 2 * z[2] / 3;
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[1][0]; p16[9][1] = p[1][1]; p16[9][2] = p[1][2];
						p16[10][0] = p[2][0]; p16[10][1] = p[2][1]; p16[10][2] = p[2][2];
						p16[11][0] = p[3][0]; p16[11][1] = p[3][1]; p16[11][2] = p[3][2];
						p16[12][0] = p[16][0]; p16[12][1] = p[16][1]; p16[12][2] = p[16][2];
						p16[13][0] = p[17][0]; p16[13][1] = p[17][1]; p16[13][2] = p[17][2];
						p16[14][0] = p[18][0]; p16[14][1] = p[18][1]; p16[14][2] = p[18][2];
						p16[15][0] = p[19][0]; p16[15][1] = p[19][1]; p16[15][2] = p[19][2];
						if (elementValenceNumber[elementValence[elementValence[i][pSId[j][1]][0]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][1]][1]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][1]][2]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[i][pSId[j][1]][3]][j][0]][5 - j] == 4) {
							p16[0][0] = 2 * p[18][0] - p[28][0] - 4 * z[0] / 3; p16[0][1] = 2 * p[18][1] - p[28][1] - 4 * z[1] / 3; p16[0][2] = 2 * p[18][2] - p[28][2] - 4 * z[2] / 3;
							p16[6][0] = 2 * p[19][0] - p[29][0] - 4 * z[0] / 3; p16[6][1] = 2 * p[19][1] - p[29][1] - 4 * z[1] / 3; p16[6][2] = 2 * p[19][2] - p[29][2] - 4 * z[2] / 3;
							p16[2][0] = p16[3][0] + p[2][0] - p[17][0]; p16[2][1] = p16[3][1] + p[2][1] - p[17][1]; p16[2][2] = p16[3][2] + p[2][2] - p[17][2];
							p16[4][0] = p16[3][0] + p[2][0] - p[16][0]; p16[4][1] = p16[3][1] + p[2][1] - p[16][1]; p16[4][2] = p16[3][2] + p[2][2] - p[16][2];
						}
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 3; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t32Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t32Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t32Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t32Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						checkOverlap = true;
						ptmp[0] = 0.5 * (p16[3][0] + p[2][0]); ptmp[1] = 0.5 * (p16[3][1] + p[2][1]); ptmp[2] = 0.5 * (p16[3][2] + p[2][2]);
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								checkOverlap = false;
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
						if (checkOverlap)// the template is already created
							hexMesh.eNum -= 3;
					}

					if (elementValenceNumber[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]] == 4) {// 3 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 15/29 6n/7n
						// 11/27 4n/5n
						// 7/23  2n/3n
						// 3/19  0n/1n
						// there are 3 hexahedrals in total
						// 0:	3	0n	2n	7	19	1n	3n	23
						// 1:	7	2n	4n	11	23	3n	5n	27
						// 2:	11	4n	6n	15	27	5n	7n	29
						// old points: 3, 7, 11, 15, 19, 23, 27, 29 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	8n	0n	2n	9n	12n	1n	3n	13n
						// 1:	9n	2n	4n	10n	13n	3n	5n	14n
						// 2:	10n	4n	6n	11n	14n	5n	7n	15n
						// points
						p16[0][0] = 2 * p[3][0] - p[2][0]; p16[0][1] = 2 * p[3][1] - p[2][1]; p16[0][2] = 2 * p[3][2] - p[2][2];
						p16[1][0] = p16[0][0] + 2 * z[0] / 3; p16[1][1] = p16[0][1] + 2 * z[1] / 3; p16[1][2] = p16[0][2] + 2 * z[2] / 3;
						p16[2][0] = 2 * p[7][0] - p[6][0]; p16[2][1] = 2 * p[7][1] - p[6][1]; p16[2][2] = 2 * p[7][2] - p[6][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[11][0] - p[10][0]; p16[4][1] = 2 * p[11][1] - p[10][1]; p16[4][2] = 2 * p[11][2] - p[10][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[6][0] = 2 * p[15][0] - p[14][0]; p16[6][1] = 2 * p[15][1] - p[14][1]; p16[6][2] = 2 * p[15][2] - p[14][2];
						p16[7][0] = p16[6][0] + 2 * z[0] / 3; p16[7][1] = p16[6][1] + 2 * z[1] / 3; p16[7][2] = p16[6][2] + 2 * z[2] / 3;
						p16[8][0] = p[3][0]; p16[8][1] = p[3][1]; p16[8][2] = p[3][2];
						p16[9][0] = p[7][0]; p16[9][1] = p[7][1]; p16[9][2] = p[7][2];
						p16[10][0] = p[11][0]; p16[10][1] = p[11][1]; p16[10][2] = p[11][2];
						p16[11][0] = p[15][0]; p16[11][1] = p[15][1]; p16[11][2] = p[15][2];
						p16[12][0] = p[19][0]; p16[12][1] = p[19][1]; p16[12][2] = p[19][2];
						p16[13][0] = p[23][0]; p16[13][1] = p[23][1]; p16[13][2] = p[23][2];
						p16[14][0] = p[27][0]; p16[14][1] = p[27][1]; p16[14][2] = p[27][2];
						p16[15][0] = p[29][0]; p16[15][1] = p[29][1]; p16[15][2] = p[29][2];
						if (elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]][0]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]][1]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]][2]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]][3]][j][0]][5 - j] == 4) {
							p16[0][0] = 2 * p[19][0] - p[18][0] - 4 * z[0] / 3; p16[0][1] = 2 * p[19][1] - p[18][1] - 4 * z[1] / 3; p16[0][2] = 2 * p[19][2] - p[18][2] - 4 * z[2] / 3;
							p16[6][0] = 2 * p[29][0] - p[28][0] - 4 * z[0] / 3; p16[6][1] = 2 * p[29][1] - p[28][1] - 4 * z[1] / 3; p16[6][2] = 2 * p[29][2] - p[28][2] - 4 * z[2] / 3;
							p16[2][0] = p16[3][0] + p[11][0] - p[27][0]; p16[2][1] = p16[3][1] + p[11][1] - p[27][1]; p16[2][2] = p16[3][2] + p[11][2] - p[27][2];
							p16[4][0] = p16[3][0] + p[11][0] - p[23][0]; p16[4][1] = p16[3][1] + p[11][1] - p[23][1]; p16[4][2] = p16[3][2] + p[11][2] - p[23][2];
						}
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 3; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t33Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t33Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t33Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t33Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						checkOverlap = true;
						ptmp[0] = 0.5 * (p16[3][0] + p[11][0]); ptmp[1] = 0.5 * (p16[3][1] + p[11][1]); ptmp[2] = 0.5 * (p16[3][2] + p[11][2]);
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								checkOverlap = false;
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
						if (checkOverlap)// the template is already created
							hexMesh.eNum -= 3;
					}

					if (elementValenceNumber[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]] == 4) {// 3 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 0n/1n 2n/3n 4n/5n 6n/7n
						// 12/28 13/30 14/31 15/29
						// there are 3 hexahedrals in total
						// 0:	12	13	2n	0n	28	30	3n	1n
						// 1:	13	14	4n	2n	30	31	5n	3n
						// 2:	14	15	6n	4n	31	29	7n	5n
						// old points: 12, 13, 14, 15, 28, 29, 30, 31 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	8n	9n	2n	0n	12n	14n	3n	1n
						// 1:	9n	10n	4n	2n	14n	15n	5n	3n
						// 2:	10n	11n	6n	4n	15n	13n	7n	5n
						// points
						p16[0][0] = 2 * p[12][0] - p[8][0]; p16[0][1] = 2 * p[12][1] - p[8][1]; p16[0][2] = 2 * p[12][2] - p[8][2];
						p16[1][0] = p16[0][0] + 2 * z[0] / 3; p16[1][1] = p16[0][1] + 2 * z[1] / 3; p16[1][2] = p16[0][2] + 2 * z[2] / 3;
						p16[2][0] = 2 * p[13][0] - p[9][0]; p16[2][1] = 2 * p[13][1] - p[9][1]; p16[2][2] = 2 * p[13][2] - p[9][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[14][0] - p[10][0]; p16[4][1] = 2 * p[14][1] - p[10][1]; p16[4][2] = 2 * p[14][2] - p[10][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[6][0] = 2 * p[15][0] - p[11][0]; p16[6][1] = 2 * p[15][1] - p[11][1]; p16[6][2] = 2 * p[15][2] - p[11][2];
						p16[7][0] = p16[6][0] + 2 * z[0] / 3; p16[7][1] = p16[6][1] + 2 * z[1] / 3; p16[7][2] = p16[6][2] + 2 * z[2] / 3;
						p16[8][0] = p[12][0]; p16[8][1] = p[12][1]; p16[8][2] = p[12][2];
						p16[9][0] = p[13][0]; p16[9][1] = p[13][1]; p16[9][2] = p[13][2];
						p16[10][0] = p[14][0]; p16[10][1] = p[14][1]; p16[10][2] = p[14][2];
						p16[11][0] = p[15][0]; p16[11][1] = p[15][1]; p16[11][2] = p[15][2];
						p16[12][0] = p[28][0]; p16[12][1] = p[28][1]; p16[12][2] = p[28][2];
						p16[13][0] = p[29][0]; p16[13][1] = p[29][1]; p16[13][2] = p[29][2];
						p16[14][0] = p[30][0]; p16[14][1] = p[30][1]; p16[14][2] = p[30][2];
						p16[15][0] = p[31][0]; p16[15][1] = p[31][1]; p16[15][2] = p[31][2];
						if (elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]][0]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]][1]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]][2]][j][0]][5 - j] == 4 ||
							elementValenceNumber[elementValence[elementValence[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]][3]][j][0]][5 - j] == 4) {
							p16[0][0] = 2 * p[28][0] - p[18][0] - 4 * z[0] / 3; p16[0][1] = 2 * p[28][1] - p[18][1] - 4 * z[1] / 3; p16[0][2] = 2 * p[28][2] - p[18][2] - 4 * z[2] / 3;
							p16[6][0] = 2 * p[29][0] - p[19][0] - 4 * z[0] / 3; p16[6][1] = 2 * p[29][1] - p[19][1] - 4 * z[1] / 3; p16[6][2] = 2 * p[29][2] - p[19][2] - 4 * z[2] / 3;
							p16[2][0] = p16[5][0] + p[13][0] - p[31][0]; p16[2][1] = p16[5][1] + p[13][1] - p[31][1]; p16[2][2] = p16[5][2] + p[13][2] - p[31][2];
							p16[4][0] = p16[5][0] + p[13][0] - p[30][0]; p16[4][1] = p16[5][1] + p[13][1] - p[30][1]; p16[4][2] = p16[5][2] + p[13][2] - p[30][2];
						}
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 3; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t34Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t34Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t34Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t34Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}

						// delete corresponding point in collectNum[][8]
						checkOverlap = true;
						ptmp[0] = 0.5 * (p16[5][0] + p[13][0]); ptmp[1] = 0.5 * (p16[5][1] + p[13][1]); ptmp[2] = 0.5 * (p16[5][2] + p[13][2]);
						for (k = 0; k < collectNumLength; k++)
							if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
								checkOverlap = false;
								collectNum.erase(collectNum.begin() + k);
								collectNumLength--;
								break;
							}
						if (checkOverlap)// the template is already created
							hexMesh.eNum -= 3;
					}

					if (elementValenceNumber[i][pSId[j][0]] == 1 && elementValenceNumber[elementValence[i][pSId[j][0]][0]][j] == 1) {// 5 element transition template

						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 6n/7n            12/28
						//       4n/5n 8/24
						//       2n/3n 4/20
						// 0n/1n            0/18
						// there are 5 hexahedrals in total
						// 0:	0n	0	4	2n	1n	18	20	3n
						// 1:	2n	4	8	4n	3n	20	24	5n
						// 2:	4n	8	12	6n	5n	24	28	7n
						// 3:	0n	2n	4n	6n	1n	3n	5n	7n
						// 4:	3n	20	24	5n	1n	18	28	7n
						// old points: 0, 4, 8, 12, 18, 20, 24, 28 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	8n	9n	2n	1n	12n	13n	3n
						// 1:	2n	9n	10n	4n	3n	13n	14n	5n
						// 2:	4n	10n	11n	6n	5n	14n	15n	7n
						// 3:	0n	2n	4n	6n	1n	3n	5n	7n
						// 4:	3n	13n	14n	5n	1n	12n	15n	7n
						// points
						p16[1][0] = 2 * p[18][0] - p[19][0]; p16[1][1] = 2 * p[18][1] - p[19][1]; p16[1][2] = 2 * p[18][2] - p[19][2];
						p16[0][0] = p16[1][0] - 4 * z[0] / 3; p16[0][1] = p16[1][1] - 4 * z[1] / 3; p16[0][2] = p16[1][2] - 4 * z[2] / 3;
						p16[2][0] = 2 * p[4][0] - p[5][0]; p16[2][1] = 2 * p[4][1] - p[5][1]; p16[2][2] = 2 * p[4][2] - p[5][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[8][0] - p[9][0]; p16[4][1] = 2 * p[8][1] - p[9][1]; p16[4][2] = 2 * p[8][2] - p[9][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[7][0] = 2 * p[28][0] - p[29][0]; p16[7][1] = 2 * p[28][1] - p[29][1]; p16[7][2] = 2 * p[28][2] - p[29][2];
						p16[6][0] = p16[7][0] - 4 * z[0] / 3; p16[6][1] = p16[7][1] - 4 * z[1] / 3; p16[6][2] = p16[7][2] - 4 * z[2] / 3;
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[4][0]; p16[9][1] = p[4][1]; p16[9][2] = p[4][2];
						p16[10][0] = p[8][0]; p16[10][1] = p[8][1]; p16[10][2] = p[8][2];
						p16[11][0] = p[12][0]; p16[11][1] = p[12][1]; p16[11][2] = p[12][2];
						p16[12][0] = p[18][0]; p16[12][1] = p[18][1]; p16[12][2] = p[18][2];
						p16[13][0] = p[20][0]; p16[13][1] = p[20][1]; p16[13][2] = p[20][2];
						p16[14][0] = p[24][0]; p16[14][1] = p[24][1]; p16[14][2] = p[24][2];
						p16[15][0] = p[28][0]; p16[15][1] = p[28][1]; p16[15][2] = p[28][2];

						p16[2][0] += p[8][0] - p[24][0] + 2 * z[0] / 3; p16[2][1] += p[8][1] - p[24][1] + 2 * z[1] / 3; p16[2][2] += p[8][2] - p[24][2] + 2 * z[2] / 3;
						p16[4][0] += p[8][0] - p[24][0] + 2 * z[0] / 3; p16[4][1] += p[8][1] - p[24][1] + 2 * z[1] / 3; p16[4][2] += p[8][2] - p[24][2] + 2 * z[2] / 3;
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 5; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t4Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t4Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t4Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t4Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}
						if (vNum == hexMesh.vNum)// overlapped element
							hexMesh.eNum -= 5;
						else {
							// delete corresponding point in collectNum[][8]
							ptmp[0] = 0.5 * (p16[2][0] + p[24][0]); ptmp[1] = 0.5 * (p16[2][1] + p[24][1]); ptmp[2] = 0.5 * (p16[2][2] + p[24][2]);
							for (k = 0; k < collectNumLength; k++)
								if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
									collectNum.erase(collectNum.begin() + k);
									collectNumLength--;
									break;
								}
						}

					}

					if (elementValenceNumber[i][pSId[j][1]] == 1 && elementValenceNumber[elementValence[i][pSId[j][1]][0]][j] == 1) {// 5 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 0/18              3/19
						//       1/16  2/17
						//       2n/3n 4n/5n
						// 0n/1n             6n/7n
						// there are 5 hexahedrals in total
						// 0:	0n	2n	1	0	1n	3n	16	18
						// 1:	2n	4n	2	1	3n	5n	17	16
						// 2:	4n	6n	3	2	5n	7n	19	17
						// 3:	3n	5n	17	16	1n	7n	19	18
						// 4:	0n	6n	4n	2n	1n	7n	5n	3n
						// old points: 0, 1, 2, 3, 16, 17, 18, 19 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	2n	9n	8n	1n	3n	12n	14n
						// 1:	2n	4n	10n	9n	3n	5n	13n	12n
						// 2:	4n	6n	11n	10n	5n	7n	15n	13n
						// 3:	3n	5n	13n	12n	1n	7n	15n	14n
						// 4:	0n	6n	4n	2n	1n	7n	5n	3n
						// points
						p16[1][0] = 2 * p[18][0] - p[28][0]; p16[1][1] = 2 * p[18][1] - p[28][1]; p16[1][2] = 2 * p[18][2] - p[28][2];
						p16[0][0] = p16[1][0] - 4 * z[0] / 3; p16[0][1] = p16[1][1] - 4 * z[1] / 3; p16[0][2] = p16[1][2] - 4 * z[2] / 3;
						p16[2][0] = 2 * p[1][0] - p[5][0]; p16[2][1] = 2 * p[1][1] - p[5][1]; p16[2][2] = 2 * p[1][2] - p[5][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[2][0] - p[6][0]; p16[4][1] = 2 * p[2][1] - p[6][1]; p16[4][2] = 2 * p[2][2] - p[6][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[7][0] = 2 * p[19][0] - p[29][0]; p16[7][1] = 2 * p[19][1] - p[29][1]; p16[7][2] = 2 * p[19][2] - p[29][2];
						p16[6][0] = p16[7][0] - 4 * z[0] / 3; p16[6][1] = p16[7][1] - 4 * z[1] / 3; p16[6][2] = p16[7][2] - 4 * z[2] / 3;
						p16[8][0] = p[0][0]; p16[8][1] = p[0][1]; p16[8][2] = p[0][2];
						p16[9][0] = p[1][0]; p16[9][1] = p[1][1]; p16[9][2] = p[1][2];
						p16[10][0] = p[2][0]; p16[10][1] = p[2][1]; p16[10][2] = p[2][2];
						p16[11][0] = p[3][0]; p16[11][1] = p[3][1]; p16[11][2] = p[3][2];
						p16[12][0] = p[16][0]; p16[12][1] = p[16][1]; p16[12][2] = p[16][2];
						p16[13][0] = p[17][0]; p16[13][1] = p[17][1]; p16[13][2] = p[17][2];
						p16[14][0] = p[18][0]; p16[14][1] = p[18][1]; p16[14][2] = p[18][2];
						p16[15][0] = p[19][0]; p16[15][1] = p[19][1]; p16[15][2] = p[19][2];

						p16[2][0] += p[1][0] - p[16][0] + 2 * z[0] / 3; p16[2][1] += p[1][1] - p[16][1] + 2 * z[1] / 3; p16[2][2] += p[1][2] - p[16][2] + 2 * z[2] / 3;
						p16[4][0] += p[1][0] - p[16][0] + 2 * z[0] / 3; p16[4][1] += p[1][1] - p[16][1] + 2 * z[1] / 3; p16[4][2] += p[1][2] - p[16][2] + 2 * z[2] / 3;
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 5; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t42Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t42Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t42Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t42Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}
						if (vNum == hexMesh.vNum)// overlapped element
							hexMesh.eNum -= 5;
						else {
							// delete corresponding point in collectNum[][8]
							ptmp[0] = 0.5 * (p16[2][0] + p[17][0]); ptmp[1] = 0.5 * (p16[2][1] + p[17][1]); ptmp[2] = 0.5 * (p16[2][2] + p[17][2]);
							for (k = 0; k < collectNumLength; k++)
								if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
									collectNum.erase(collectNum.begin() + k);
									collectNumLength--;
									break;
								}
						}
					}

					if (elementValenceNumber[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]] == 1 &&
						elementValenceNumber[elementValence[elementValence[i][pS2Id[j][0]][0]][pS2Id[j][0]][0]][j] == 1) {// 5 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 15/29             6n/7n
						//       11/27 4n/5n 
						//       7/23  2n/3n
						// 3/19              0n/1n
						// there are 5 hexahedrals in total
						// 0:	3	0n	2n	7	19	1n	3n	23
						// 1:	7	2n	4n	11	23	3n	5n	27
						// 2:	11	4n	6n	15	27	5n	7n	29
						// 3:	23	3n	5n	27	19	1n	7n	29
						// 4:	2n	0n	6n	4n	3n	1n	7n	5n
						// old points: 3, 7, 11, 15, 19, 23, 27, 29 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	8n	0n	2n	9n	12n	1n	3n	13n
						// 1:	9n	2n	4n	10n	13n	3n	5n	14n
						// 2:	10n	4n	6n	11n	14n	5n	7n	15n
						// 3:	13n	3n	5n	14n	12n	1n	7n	15n
						// 4:	2n	0n	6n	4n	3n	1n	7n	5n
						// points
						p16[1][0] = 2 * p[19][0] - p[18][0]; p16[1][1] = 2 * p[19][1] - p[18][1]; p16[1][2] = 2 * p[19][2] - p[18][2];
						p16[0][0] = p16[1][0] - 4 * z[0] / 3; p16[0][1] = p16[1][1] - 4 * z[1] / 3; p16[0][2] = p16[1][2] - 4 * z[2] / 3;
						p16[2][0] = 2 * p[7][0] - p[6][0]; p16[2][1] = 2 * p[7][1] - p[6][1]; p16[2][2] = 2 * p[7][2] - p[6][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[11][0] - p[10][0]; p16[4][1] = 2 * p[11][1] - p[10][1]; p16[4][2] = 2 * p[11][2] - p[10][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[7][0] = 2 * p[29][0] - p[28][0]; p16[7][1] = 2 * p[29][1] - p[28][1]; p16[7][2] = 2 * p[29][2] - p[28][2];
						p16[6][0] = p16[7][0] - 4 * z[0] / 3; p16[6][1] = p16[7][1] - 4 * z[1] / 3; p16[6][2] = p16[7][2] - 4 * z[2] / 3;
						p16[8][0] = p[3][0]; p16[8][1] = p[3][1]; p16[8][2] = p[3][2];
						p16[9][0] = p[7][0]; p16[9][1] = p[7][1]; p16[9][2] = p[7][2];
						p16[10][0] = p[11][0]; p16[10][1] = p[11][1]; p16[10][2] = p[11][2];
						p16[11][0] = p[15][0]; p16[11][1] = p[15][1]; p16[11][2] = p[15][2];
						p16[12][0] = p[19][0]; p16[12][1] = p[19][1]; p16[12][2] = p[19][2];
						p16[13][0] = p[23][0]; p16[13][1] = p[23][1]; p16[13][2] = p[23][2];
						p16[14][0] = p[27][0]; p16[14][1] = p[27][1]; p16[14][2] = p[27][2];
						p16[15][0] = p[29][0]; p16[15][1] = p[29][1]; p16[15][2] = p[29][2];

						p16[2][0] += p[7][0] - p[23][0] + 2 * z[0] / 3; p16[2][1] += p[7][1] - p[23][1] + 2 * z[1] / 3; p16[2][2] += p[7][2] - p[23][2] + 2 * z[2] / 3;
						p16[4][0] += p[7][0] - p[23][0] + 2 * z[0] / 3; p16[4][1] += p[7][1] - p[23][1] + 2 * z[1] / 3; p16[4][2] += p[7][2] - p[23][2] + 2 * z[2] / 3;
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 5; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t43Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t43Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t43Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t43Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}
						if (vNum == hexMesh.vNum)// overlapped element
							hexMesh.eNum -= 5;
						else {
							// delete corresponding point in collectNum[][8]
							ptmp[0] = 0.5 * (p16[2][0] + p[27][0]); ptmp[1] = 0.5 * (p16[2][1] + p[27][1]); ptmp[2] = 0.5 * (p16[2][2] + p[27][2]);
							for (k = 0; k < collectNumLength; k++)
								if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
									collectNum.erase(collectNum.begin() + k);
									collectNumLength--;
									break;
								}
						}
					}

					if (elementValenceNumber[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]] == 1 &&
						elementValenceNumber[elementValence[elementValence[i][pS2Id[j][1]][0]][pS2Id[j][1]][0]][j] == 1) {// 5 element transition template
						// there are 16 vertices in total, 8 of them are new (marked as *n below)
						// 0n/1n             6n/7n
						//       2n/3n 4n/5n 
						//       13/30 14/31
						// 12/28             15/29
						// there are 5 hexahedrals in total
						// 0:	0n	12	13	2n	1n	28	30	3n
						// 1:	2n	13	14	4n	3n	30	31	5n
						// 2:	4n	14	15	6n	5n	31	29	7n
						// 3:	3n	30	31	5n	1n	28	29	7n
						// 4:	0n	2n	4n	6n	1n	3n	5n	7n
						// old points: 12, 13, 14, 15, 28, 29, 30, 31 -> new points: 8n, 9n, 10n, 11n, 12n, 13n, 14n, 15n
						// 0:	0n	8n	9n	2n	1n	12n	14n	3n
						// 1:	2n	9n	10n	4n	3n	14n	15n	5n
						// 2:	4n	10n	11n	6n	5n	15n	13n	7n
						// 3:	3n	14n	15n	5n	1n	12n	13n	7n
						// 4:	0n	2n	4n	6n	1n	3n	5n	7n
						// points
						p16[1][0] = 2 * p[28][0] - p[18][0]; p16[1][1] = 2 * p[28][1] - p[18][1]; p16[1][2] = 2 * p[28][2] - p[18][2];
						p16[0][0] = p16[1][0] - 4 * z[0] / 3; p16[0][1] = p16[1][1] - 4 * z[1] / 3; p16[0][2] = p16[1][2] - 4 * z[2] / 3;
						p16[2][0] = 2 * p[13][0] - p[9][0]; p16[2][1] = 2 * p[13][1] - p[9][1]; p16[2][2] = 2 * p[13][2] - p[9][2];
						p16[3][0] = p16[2][0] + 2 * z[0] / 3; p16[3][1] = p16[2][1] + 2 * z[1] / 3; p16[3][2] = p16[2][2] + 2 * z[2] / 3;
						p16[4][0] = 2 * p[14][0] - p[10][0]; p16[4][1] = 2 * p[14][1] - p[10][1]; p16[4][2] = 2 * p[14][2] - p[10][2];
						p16[5][0] = p16[4][0] + 2 * z[0] / 3; p16[5][1] = p16[4][1] + 2 * z[1] / 3; p16[5][2] = p16[4][2] + 2 * z[2] / 3;
						p16[7][0] = 2 * p[29][0] - p[19][0]; p16[7][1] = 2 * p[29][1] - p[19][1]; p16[7][2] = 2 * p[29][2] - p[19][2];
						p16[6][0] = p16[7][0] - 4 * z[0] / 3; p16[6][1] = p16[7][1] - 4 * z[1] / 3; p16[6][2] = p16[7][2] - 4 * z[2] / 3;
						p16[8][0] = p[12][0]; p16[8][1] = p[12][1]; p16[8][2] = p[12][2];
						p16[9][0] = p[13][0]; p16[9][1] = p[13][1]; p16[9][2] = p[13][2];
						p16[10][0] = p[14][0]; p16[10][1] = p[14][1]; p16[10][2] = p[14][2];
						p16[11][0] = p[15][0]; p16[11][1] = p[15][1]; p16[11][2] = p[15][2];
						p16[12][0] = p[28][0]; p16[12][1] = p[28][1]; p16[12][2] = p[28][2];
						p16[13][0] = p[29][0]; p16[13][1] = p[29][1]; p16[13][2] = p[29][2];
						p16[14][0] = p[30][0]; p16[14][1] = p[30][1]; p16[14][2] = p[30][2];
						p16[15][0] = p[31][0]; p16[15][1] = p[31][1]; p16[15][2] = p[31][2];

						p16[2][0] += p[13][0] - p[30][0] + 2 * z[0] / 3; p16[2][1] += p[13][1] - p[30][1] + 2 * z[1] / 3; p16[2][2] += p[13][2] - p[30][2] + 2 * z[2] / 3;
						p16[4][0] += p[13][0] - p[30][0] + 2 * z[0] / 3; p16[4][1] += p[13][1] - p[30][1] + 2 * z[1] / 3; p16[4][2] += p[13][2] - p[30][2] + 2 * z[2] / 3;
						// elements
						vNum = hexMesh.vNum;
						for (k = 0; k < 5; k++) {
							for (l = 0; l < 8; l++) {
								overlap = false;
								for (m = 0; m < hexMesh.vNum; m++)
									if (dist(hexMesh.v[m], p16[t44Id[stepI][k][l]]) < DIST_THRES) {// overlap
										overlap = true;
										hexMesh.e[hexMesh.eNum][l] = m;
										break;
									}
								if (!overlap) {// create new point
									hexMesh.e[hexMesh.eNum][l] = hexMesh.vNum;
									hexMesh.v[hexMesh.vNum][0] = p16[t44Id[stepI][k][l]][0];
									hexMesh.v[hexMesh.vNum][1] = p16[t44Id[stepI][k][l]][1];
									hexMesh.v[hexMesh.vNum][2] = p16[t44Id[stepI][k][l]][2];
									hexMesh.vNum++;
								}
							}
							hexMesh.eNum++;
						}
						if (vNum == hexMesh.vNum)// overlapped element
							hexMesh.eNum -= 5;
						else {
							// delete corresponding point in collectNum[][8]
							ptmp[0] = 0.5 * (p16[2][0] + p[31][0]); ptmp[1] = 0.5 * (p16[2][1] + p[31][1]); ptmp[2] = 0.5 * (p16[2][2] + p[31][2]);
							for (k = 0; k < collectNumLength; k++)
								if (dist(octreeMesh.v[collectNum[k][8]], ptmp) < DIST_THRES) {
									collectNum.erase(collectNum.begin() + k);
									collectNumLength--;
									break;
								}
						}
					}
				}
			}
	}

	for (i = 0; i < collectNumLength; i++) {
		for (j = 0; j < 8; j++) {
			overlap = false;
			for (k = 0; k < hexMesh.vNum; k++)
				if (dist(hexMesh.v[k], tmp[collectNum[i][j]]) < DIST_THRES) {// overlap
					overlap = true;
					hexMesh.e[hexMesh.eNum][j] = k;
					break;
				}
			if (!overlap) {// create new point
				hexMesh.e[hexMesh.eNum][j] = hexMesh.vNum;
				hexMesh.v[hexMesh.vNum][0] = tmp[collectNum[i][j]][0];
				hexMesh.v[hexMesh.vNum][1] = tmp[collectNum[i][j]][1];
				hexMesh.v[hexMesh.vNum][2] = tmp[collectNum[i][j]][2];
				hexMesh.vNum++;
			}
		}
		hexMesh.eNum++;
	}

	for (i = 0; i < hexMesh.vNum; i++) {
		hexMesh.v[i][0] = START_POINT[0] + BOX_LENGTH_RATIO * hexMesh.v[i][0];
		hexMesh.v[i][1] = START_POINT[1] + BOX_LENGTH_RATIO * hexMesh.v[i][1];
		hexMesh.v[i][2] = START_POINT[2] + BOX_LENGTH_RATIO * hexMesh.v[i][2];
	}

	BOX_LENGTH_RATIO = 1;  START_POINT[0] = 0; START_POINT[1] = 0; START_POINT[2] = 0;

	hexMesh.WriteToVtk(fileName, BOX_LENGTH_RATIO, START_POINT);
}

void hexGen::ReadDualFullHex(const char* inputFileName) {
	BOX_LENGTH_RATIO = 1;  START_POINT[0] = 0; START_POINT[1] = 0; START_POINT[2] = 0;
	FILE* dataFile = fopen(inputFileName, "r");
	if (NULL == dataFile)
	{
		std::cerr << "ErrorCode 0: Wrong file name " << inputFileName << std::endl;
		return;
	}

	char line[256];
	int i, j;
	for (i = 0; i < 4 && fgets(line, sizeof(line), dataFile); i++) {
	}

	int points, elements;
	if (fgets(line, sizeof(line), dataFile) && sscanf(line, "POINTS %d double", &points) == 1) {
		hexMesh.Initialize(points);
		hexMeshENum = points;
		hexMesh.vNum = points;
		for (i = 0; i < points; i++) {
			fgets(line, sizeof(line), dataFile);
			sscanf(line, "%lf %lf %lf", &hexMesh.v[i][0], &hexMesh.v[i][1], &hexMesh.v[i][2]);
		}
		if (fgets(line, sizeof(line), dataFile) && sscanf(line, "CELLS %d %d", &elements, &i) == 2) {
			hexMesh.eNum = elements;
			for (i = 0; i < elements; i++) {
				fgets(line, sizeof(line), dataFile);
				sscanf(line, "%d %d %d %d %d %d %d %d %d", &j, &hexMesh.e[i][0], &hexMesh.e[i][1], &hexMesh.e[i][2], &hexMesh.e[i][3], &hexMesh.e[i][4], &hexMesh.e[i][5], &hexMesh.e[i][6], &hexMesh.e[i][7]);
			}
		}
		else
			std::cerr << "ErrorCode 2: Cannot get total element number" << std::endl;
	}
	else
		std::cerr << "ErrorCode 1: Cannot get total point number" << std::endl;
	fclose(dataFile);
}

void hexGen::RemoveOutsideElement(const char* fileName) {// write from hexMesh to octreeMesh (leafNum)
	int i, j, k, l, idx[27]; bool pass, inside, overlap; double alpha, ref[19][3], tmp[3], tmp2[3], dir[3];

	leafNum = 0; octreeMesh.vNum = 0; std::vector<double> deletePoint(hexMesh.vNum, MAX_NUM2);

	for (i = 0; i < hexMesh.vNum; i++) {
		pass = false;
		while (!pass) {
			pass = true;  inside = false;
			tmp2[0] = (rand() / (RAND_MAX * 1.f) + DIST_THRES) * ((rand() / (RAND_MAX * 1.f) - 0.5 > 0) ? -1 : 1);
			tmp2[1] = (rand() / (RAND_MAX * 1.f) + DIST_THRES) * ((rand() / (RAND_MAX * 1.f) - 0.5 > 0) ? -1 : 1);
			tmp2[2] = (rand() / (RAND_MAX * 1.f) + DIST_THRES) * ((rand() / (RAND_MAX * 1.f) - 0.5 > 0) ? -1 : 1);
			for (j = 0; j < triMesh.eNum; j++) {
				k = Intersect(triMesh.v[triMesh.e[j][0]], triMesh.v[triMesh.e[j][1]], triMesh.v[triMesh.e[j][2]], hexMesh.v[i], tmp2, tmp, alpha);
				if (k == 1 && alpha > 0) inside = !inside;
				else if (k == -1) {
					pass = false; break;
				}
			}
		}
		for (j = 0; j < triMesh.eNum; j++) {
			alpha = PointToTri(triMesh.v[triMesh.e[j][0]], triMesh.v[triMesh.e[j][1]], triMesh.v[triMesh.e[j][2]], hexMesh.v[i], tmp, deletePoint[i]);
			if (alpha < deletePoint[i])
				deletePoint[i] = alpha;
		}
		if (!inside) deletePoint[i] = -deletePoint[i];// signed distance, negative for outside points
	}

	for (i = 0; i < hexMesh.eNum; i++) {
		tmp[0] = 0; tmp[1] = 0; k = 0;
		for (j = 0; j < 8; j++) {
			if (deletePoint[hexMesh.e[i][j]] > tmp[0]) tmp[0] = deletePoint[hexMesh.e[i][j]];
			else if (deletePoint[hexMesh.e[i][j]] < 0) {
				k++; if (k > 2) break;
				if (deletePoint[hexMesh.e[i][j]] < tmp[1]) tmp[1] = deletePoint[hexMesh.e[i][j]];
			}
		}
		if (k < 3 && tmp[1] + OUT_IN_RATIO * tmp[0] >= 0) {
			for (j = 0; j < 8; j++) {
				overlap = false;
				for (k = 0; k < octreeMesh.vNum; k++)
					if (dist(octreeMesh.v[k], hexMesh.v[hexMesh.e[i][j]]) < DIST_THRES) {// overlap
						overlap = true;
						octreeMesh.e[leafNum][j] = k;
						break;
					}
				if (!overlap) {// create new point
					octreeMesh.e[leafNum][j] = octreeMesh.vNum;
					octreeMesh.v[octreeMesh.vNum][0] = hexMesh.v[hexMesh.e[i][j]][0];
					octreeMesh.v[octreeMesh.vNum][1] = hexMesh.v[hexMesh.e[i][j]][1];
					octreeMesh.v[octreeMesh.vNum][2] = hexMesh.v[hexMesh.e[i][j]][2];
					octreeMesh.vNum++;
				}
			}
			leafNum++;
		}
	}
	octreeMesh.eNum = leafNum; octreeMesh.WriteToVtk(fileName, BOX_LENGTH_RATIO, START_POINT);

	// delete invalid topologies
	overlap = true; int fIdx = 0;
	std::vector<bool> deleteElement(leafNum, false);
	int minIdx, minF, m, m0, m1, m2, n0, n1, n2, faceTmp[6][5], faceTmp2[6][5];
	double x1[3], x2[3];
	std::vector <int> toBeDelete;

	std::vector<std::vector<bool>> ck(octreeMesh.vNum, std::vector<bool>(146, true));
	std::vector<std::vector<int>> face(leafNum * 6, std::vector<int>(5));// 4 points in surface + exist number
	// initialize face overlap number
	for (i = 0; i < leafNum * 6; i++) face[i][4] = -1;// -1: unused; [0, leafNum): element number

	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 6; j++) {
			inside = false;
			for (l = 0; l < fIdx; l++) {// put k inside for more convenient break operation
				for (k = 0; k < 4; k++)
					if (face[l][0] == octreeMesh.e[i][fIdC[j][k]] && face[l][2] == octreeMesh.e[i][fIdC[j][(k + 2) % 4]]) {// find two matching faces
						inside = true; fIdx--; face.erase(face.begin() + l);
						break;
					}
				if (inside) break;
			}
			if (!inside) {
				face[fIdx][4] = i;
				face[fIdx][0] = octreeMesh.e[i][fIdC[j][0]]; face[fIdx][1] = octreeMesh.e[i][fIdC[j][1]]; face[fIdx][2] = octreeMesh.e[i][fIdC[j][2]]; face[fIdx][3] = octreeMesh.e[i][fIdC[j][3]];
				fIdx++;
			}
		}

	while (true) {
		for (i = 0; i < octreeMesh.vNum; i++)
			for (j = 0; j < 146; j++) ck[i][j] = true;
		for (i = 0; i < fIdx; i++)
			for (j = 0; j < 4; j++) {
				n0 = face[i][j]; n1 = face[i][(j == 3) ? 0 : j + 1]; n2 = face[i][(j == 0) ? 3 : j - 1];
				//need to judge direction is toward inside or outside
				for (k = 0; k < 8; k++) {
					if (n0 == octreeMesh.e[face[i][4]][k]) m0 = k;
					if (n1 == octreeMesh.e[face[i][4]][k]) m1 = k;
					if (n2 == octreeMesh.e[face[i][4]][k]) m2 = k;
				}
				if (!towardOutside[m0][m1][m2]) {
					m0 = n1; n1 = n2; n2 = m0;
				}

				x1[0] = octreeMesh.v[n1][0] - octreeMesh.v[n0][0];
				x1[1] = octreeMesh.v[n1][1] - octreeMesh.v[n0][1];
				x1[2] = octreeMesh.v[n1][2] - octreeMesh.v[n0][2];
				x2[0] = octreeMesh.v[n2][0] - octreeMesh.v[n0][0];
				x2[1] = octreeMesh.v[n2][1] - octreeMesh.v[n0][1];
				x2[2] = octreeMesh.v[n2][2] - octreeMesh.v[n0][2];

				tmp[0] = x1[1] * x2[2] - x1[2] * x2[1];
				tmp[1] = x1[2] * x2[0] - x1[0] * x2[2];
				tmp[2] = x1[0] * x2[1] - x1[1] * x2[0];
				for (k = 0; k < 146; k++)
					if (ck[n0][k] && DOT(tmp, pointOnSurf[k]) < COS_THRES * sqrt(dist(0, 0, 0, tmp))) ck[n0][k] = false;
			}
		minF = 6;
		for (i = 0; i < fIdx; i++)
			for (j = 0; j < 4; j++) {
				pass = false; n1 = face[i][j];
				for (k = 0; k < 146; k++)
					if (ck[n1][k]) {
						pass = true; break;
					}
				if (pass) continue;

				m = face[i][4]; n0 = 0; deleteElement[m] = true;
				for (l = 0; l < 6; l++) {
					for (m2 = 0; m2 < 4; m2++)
						faceTmp[l][m2] = octreeMesh.e[m][fId[l][m2]];
					faceTmp[l][4] = -1;
					overlap = false;
					for (m2 = 0; m2 < leafNum; m2++)
						if (!deleteElement[m2]) {
							for (m0 = 0; m0 < 6; m0++) {
								for (m1 = 0; m1 < 4; m1++)
									if (faceTmp[l][0] == octreeMesh.e[m2][fId[m0][m1]] && faceTmp[l][2] == octreeMesh.e[m2][fId[m0][(m1 + 2) % 4]]) {
										faceTmp[l][4] = m2; overlap = true; break;
									}
								if (overlap) break;
							}
							if (overlap) break;
						}
					if (overlap) n0++;
				}
				if (n0 < minF) {
					minF = n0; minIdx = i;
					for (l = 0; l < 6; l++)
						for (m2 = 0; m2 < 5; m2++) faceTmp2[l][m2] = faceTmp[l][m2];
				}
				deleteElement[m] = false;
			}
		if (minF == 6) break;
		i = minIdx; deleteElement[face[i][4]] = true;
		std::cout << face[i][4] << ", ";

		for (l = fIdx - 1; l > -1; l--)
			if (face[l][4] == face[i][4]) toBeDelete.push_back(l);
		for (int l : toBeDelete) {
			face.erase(face.begin() + l);
		}
		fIdx -= toBeDelete.size();
		toBeDelete.clear();

		for (l = 0; l < 6; l++)
			if (faceTmp2[l][4] != -1) {
				for (m0 = 0; m0 < 5; m0++) face[fIdx][m0] = faceTmp2[l][m0];
				std::cout << face[fIdx][0] << " " << face[fIdx][1] << " " << face[fIdx][2] << " " << face[fIdx][3] << " " << face[fIdx][4] << " " << std::endl;
				fIdx++;
			}
	}
	l = leafNum; leafNum = 0; hexMesh.vNum = 0;

	for (i = 0; i < l; i++) {
		if (deleteElement[i]) continue;

		for (j = 0; j < 8; j++) {
			overlap = false;
			for (k = 0; k < hexMesh.vNum; k++)
				if (dist(hexMesh.v[k], octreeMesh.v[octreeMesh.e[i][j]]) < DIST_THRES) {// overlap
					overlap = true;
					hexMesh.e[leafNum][j] = k;
					break;
				}
			if (!overlap) {// create new point
				hexMesh.e[leafNum][j] = hexMesh.vNum;
				hexMesh.v[hexMesh.vNum][0] = octreeMesh.v[octreeMesh.e[i][j]][0];
				hexMesh.v[hexMesh.vNum][1] = octreeMesh.v[octreeMesh.e[i][j]][1];
				hexMesh.v[hexMesh.vNum][2] = octreeMesh.v[octreeMesh.e[i][j]][2];
				hexMesh.vNum++;
			}
		}
		leafNum++;
	}
	octreeMesh.vNum = hexMesh.vNum;  octreeMesh.eNum = leafNum;
	for (i = 0; i < octreeMesh.vNum; i++) {
		octreeMesh.v[i][0] = hexMesh.v[i][0]; octreeMesh.v[i][1] = hexMesh.v[i][1]; octreeMesh.v[i][2] = hexMesh.v[i][2];
	}
	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 8; j++)
			octreeMesh.e[i][j] = hexMesh.e[i][j];

	octreeMesh.WriteToVtk(fileName, BOX_LENGTH_RATIO, START_POINT);
}

void hexGen::ReadDualHex(const char* inputFileName) {
	// here octreeMesh is used again, but this is not octreeMesh, it is just saving memory space
	// leafNum is always with octreeMesh, but not with hexMesh
	FILE* dataFile = fopen(inputFileName, "r");
	if (NULL == dataFile)
	{
		std::cerr << "ErrorCode 0: Wrong file name " << inputFileName << std::endl;
		return;
	}

	char line[256];
	int i, j;
	for (i = 0; i < 4 && fgets(line, sizeof(line), dataFile); i++) {
	}

	int points;
	if (fgets(line, sizeof(line), dataFile) && sscanf(line, "POINTS %d double", &points) == 1) {
		delete[] octreeMesh.v;
		octreeMesh.v = nullptr;
		for (i = 0; i < octreeENum; i++) {
			delete[] octreeMesh.e[i];
			octreeMesh.e[i] = nullptr;
		}
		delete[] octreeMesh.e;
		octreeMesh.e = nullptr;
		octreeMesh.Initialize(points * 2, points * 2);
		octreeMesh.vNum = points;
		for (i = 0; i < points; i++) {
			fgets(line, sizeof(line), dataFile);
			sscanf(line, "%lf %lf %lf", &octreeMesh.v[i][0], &octreeMesh.v[i][1], &octreeMesh.v[i][2]);
		}
		if (fgets(line, sizeof(line), dataFile) && sscanf(line, "CELLS %d %d", &leafNum, &i) == 2) {
			octreeMesh.eNum = leafNum;
			for (i = 0; i < leafNum; i++) {
				fgets(line, sizeof(line), dataFile);
				sscanf(line, "%d %d %d %d %d %d %d %d %d", &j, &octreeMesh.e[i][0], &octreeMesh.e[i][1], &octreeMesh.e[i][2], &octreeMesh.e[i][3], &octreeMesh.e[i][4], &octreeMesh.e[i][5], &octreeMesh.e[i][6], &octreeMesh.e[i][7]);
			}
		}
		else
			std::cerr << "ErrorCode 2: Cannot get total element number" << std::endl;
	}
	else
		std::cerr << "ErrorCode 1: Cannot get total point number" << std::endl;
	fclose(dataFile);
}

void hexGen::ProjectToIsoSurface(const char* fileName) {// modify octreeMesh only
	int i, j, k, l, fIdx = 0; bool pair;
	std::vector<std::vector<int>> face(leafNum * 6, std::vector<int>(5));// 4 points in surface + exist number
	// initialize face overlap number
	for (i = 0; i < leafNum * 6; i++) face[i][4] = -1;// -1: unused; [0, leafNum): element number

	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 6; j++) {
			pair = false;
			for (l = 0; l < fIdx; l++) {// put k inside for more convenient break operation
				for (k = 0; k < 4; k++)
					if (face[l][0] == octreeMesh.e[i][fIdC[j][k]] && face[l][2] == octreeMesh.e[i][fIdC[j][(k + 2) % 4]] && (
						(face[l][1] == octreeMesh.e[i][fIdC[j][(k + 1) % 4]] && face[l][3] == octreeMesh.e[i][fIdC[j][(k + 3) % 4]]) ||
						(face[l][3] == octreeMesh.e[i][fIdC[j][(k + 1) % 4]] && face[l][1] == octreeMesh.e[i][fIdC[j][(k + 3) % 4]]))) {// find two matching faces
						pair = true; fIdx--; face.erase(face.begin() + l);
						break;
					}
				if (pair) break;
			}
			if (!pair) {
				face[fIdx][4] = i;
				face[fIdx][0] = octreeMesh.e[i][fIdC[j][0]]; face[fIdx][1] = octreeMesh.e[i][fIdC[j][1]]; face[fIdx][2] = octreeMesh.e[i][fIdC[j][2]]; face[fIdx][3] = octreeMesh.e[i][fIdC[j][3]];
				fIdx++;
			}
		}
	face.erase(face.begin() + fIdx, face.end());

	std::vector<bool> pointOnSurface(octreeMesh.vNum, false);// if a point is on the surface
	std::vector<int> projNum(octreeMesh.vNum);
	int sPIdx = 0;

	for (i = 0; i < fIdx; i++)
		for (j = 0; j < 4; j++)
			pointOnSurface[face[i][j]] = true;
	for (i = 0; i < octreeMesh.vNum; i++)
		if (pointOnSurface[i]) sPIdx++;
	std::vector<int> bP(2 * sPIdx);
	for (i = 0; i < sPIdx; i++)
		bP[i + sPIdx] = i + octreeMesh.vNum;
	sPIdx = 0;
	for (i = 0; i < octreeMesh.vNum; i++)
		if (pointOnSurface[i]) {
			bP[sPIdx] = i;
			projNum[i] = sPIdx + octreeMesh.vNum;
			sPIdx++;
		}

	pointOnSurface.clear();

	std::vector<std::vector<int>> aE(leafNum + fIdx, std::vector<int>(8, -1));
	std::vector<int> bP2;
	std::vector<std::vector<int>> p2E(octreeMesh.vNum, std::vector<int>());
	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 8; j++)
			p2E[octreeMesh.e[i][j]].push_back(i);
	for (i = 0; i < sPIdx; i++)
		for (j = 0; j < p2E[bP[i]].size(); j++)
			for (k = 0; k < 8; k++)
				if (octreeMesh.e[p2E[bP[i]][j]][k] == bP[i]) {
					aE[p2E[bP[i]][j]][k] = i; break;
				}

	int affElemNum = leafNum;
	for (i = 0; i < octreeMesh.vNum; i++)
		for (j = 0; j < p2E[i].size(); j++) {
			for (k = 0; k < 8; k++)
				if (aE[p2E[i][j]][k] == -1 && octreeMesh.e[p2E[i][j]][k] == i) {
					aE[p2E[i][j]][k] = 2 * sPIdx + bP2.size();
					bP2.push_back(i);
					break;
				}
		}

	for (i = 0; i < sPIdx; i++) {
		octreeMesh.v[i + octreeMesh.vNum][0] = octreeMesh.v[bP[i]][0];
		octreeMesh.v[i + octreeMesh.vNum][1] = octreeMesh.v[bP[i]][1];
		octreeMesh.v[i + octreeMesh.vNum][2] = octreeMesh.v[bP[i]][2];
	}
	for (i = 0; i < fIdx; i++) {
		for (j = 0; j < 4; j++) {
			octreeMesh.e[leafNum][j + 4] = face[i][j];
			for (k = 0; k < sPIdx; k++)
				if (bP[k] == face[i][j]) {
					octreeMesh.e[leafNum][j] = bP[k + sPIdx];
					aE[affElemNum][j] = k + sPIdx;
					aE[affElemNum][j + 4] = k;
					break;
				}
		}
		affElemNum++;
		leafNum++;
	}
	octreeMesh.vNum += sPIdx;
	octreeMesh.eNum = leafNum;

	std::vector<std::vector<int>> cE(sPIdx);
	std::vector<std::vector<int>> cP(sPIdx);
	std::vector<std::vector<int>> cE2(sPIdx);
	std::vector<std::vector<int>> cP2(sPIdx);
	double testP[8][3];

	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 8; j++)
			for (k = 0; k < sPIdx; k++) {
				if (octreeMesh.e[i][j] == bP[k + sPIdx]) {
					cE[k].push_back(i);
					break;
				}

				if (octreeMesh.e[i][j] == bP[k]) {
					cE2[k].push_back(i);

					for (l = 0; l < cP2[k].size(); l++)
						if (cP2[k][l] == octreeMesh.e[i][adjP[j][0]]) break;

					if (l == cP2[k].size())
						cP2[k].push_back(octreeMesh.e[i][adjP[j][0]]);

					for (l = 0; l < cP2[k].size(); l++)
						if (cP2[k][l] == octreeMesh.e[i][adjP[j][1]]) break;

					if (l == cP2[k].size())
						cP2[k].push_back(octreeMesh.e[i][adjP[j][1]]);

					for (l = 0; l < cP2[k].size(); l++)
						if (cP2[k][l] == octreeMesh.e[i][adjP[j][2]]) break;

					if (l == cP2[k].size())
						cP2[k].push_back(octreeMesh.e[i][adjP[j][2]]);
					break;
				}
			}

	for (i = 0; i < fIdx; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < sPIdx; k++)
				if (face[i][j] == bP[k]) {

					for (l = 0; l < cP[k].size(); l++)
						if (cP[k][l] == projNum[face[i][(j == 0) ? 3 : (j - 1)]]) break;

					if (l == cP[k].size())
						cP[k].push_back(projNum[face[i][(j == 0) ? 3 : (j - 1)]]);

					for (l = 0; l < cP[k].size(); l++)
						if (cP[k][l] == projNum[face[i][(j == 3) ? 0 : (j + 1)]]) break;

					if (l == cP[k].size())
						cP[k].push_back(projNum[face[i][(j == 3) ? 0 : (j + 1)]]);
					break;
				}
	projNum.clear();

	double(*g)[3], target[3], minDist, dis, sJL = 0, tmp[3], tmp2[3], LEARNING_RATE = 5.0e-4f, aveDist = MAX_NUM2, smallDist = 114514, prop,
		x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7;
	g = new double[2 * sPIdx + bP2.size()][3];// gradient
	int minIdx[9], maxDistIdx = -114514;
	bool allPositive = false;

	// choose which triangles to project to for all points
	std::vector<int> triNum(sPIdx);
	// brute force drag count
	std::vector<int> dragCount(sPIdx, 0);

	for (i = 0; i < sPIdx; i++) {
		minDist = MAX_NUM2;
		for (j = 0; j < triMesh.eNum; j++) {
			dis = PointToTri(triMesh.v[triMesh.e[j][0]], triMesh.v[triMesh.e[j][1]], triMesh.v[triMesh.e[j][2]], octreeMesh.v[bP[i]], target, minDist);
			if (dis < minDist) {
				minDist = dis;
				triNum[i] = j;
			}
		}
	}
	i = 0;
	while (true) {
		i++;
		for (j = 0; j < 2 * sPIdx + bP2.size(); j++) {
			g[j][0] = 0; g[j][1] = 0; g[j][2] = 0;
		}

		// optimize scaled jacobian
		sJL = 0;
		for (j = 0; j < affElemNum; j++) {
			// calculate bad jacobian number
			x0 = octreeMesh.v[octreeMesh.e[j][0]][0];
			y0 = octreeMesh.v[octreeMesh.e[j][0]][1];
			z0 = octreeMesh.v[octreeMesh.e[j][0]][2];
			x1 = octreeMesh.v[octreeMesh.e[j][1]][0];
			y1 = octreeMesh.v[octreeMesh.e[j][1]][1];
			z1 = octreeMesh.v[octreeMesh.e[j][1]][2];
			x2 = octreeMesh.v[octreeMesh.e[j][2]][0];
			y2 = octreeMesh.v[octreeMesh.e[j][2]][1];
			z2 = octreeMesh.v[octreeMesh.e[j][2]][2];
			x3 = octreeMesh.v[octreeMesh.e[j][3]][0];
			y3 = octreeMesh.v[octreeMesh.e[j][3]][1];
			z3 = octreeMesh.v[octreeMesh.e[j][3]][2];
			x4 = octreeMesh.v[octreeMesh.e[j][4]][0];
			y4 = octreeMesh.v[octreeMesh.e[j][4]][1];
			z4 = octreeMesh.v[octreeMesh.e[j][4]][2];
			x5 = octreeMesh.v[octreeMesh.e[j][5]][0];
			y5 = octreeMesh.v[octreeMesh.e[j][5]][1];
			z5 = octreeMesh.v[octreeMesh.e[j][5]][2];
			x6 = octreeMesh.v[octreeMesh.e[j][6]][0];
			y6 = octreeMesh.v[octreeMesh.e[j][6]][1];
			z6 = octreeMesh.v[octreeMesh.e[j][6]][2];
			x7 = octreeMesh.v[octreeMesh.e[j][7]][0];
			y7 = octreeMesh.v[octreeMesh.e[j][7]][1];
			z7 = octreeMesh.v[octreeMesh.e[j][7]][2];
			iSj(octreeMesh.v[octreeMesh.e[j][0]], octreeMesh.v[octreeMesh.e[j][1]],
				octreeMesh.v[octreeMesh.e[j][2]], octreeMesh.v[octreeMesh.e[j][3]],
				octreeMesh.v[octreeMesh.e[j][4]], octreeMesh.v[octreeMesh.e[j][5]],
				octreeMesh.v[octreeMesh.e[j][6]], octreeMesh.v[octreeMesh.e[j][7]], minIdx);
			pair = false;
			if (minIdx[0] == 0) {
				pair = true;
				g[aE[j][0]][0] += -4 * (y4 * z1 + y5 * z1 + y1 * z2 + y5 * z2 - y7 * z2 + y1 * z3 - y4 * z3 - y7 * z3 -
					y1 * z4 - y5 * z4 + y7 * z4 - y1 * z5 + y4 * z5 + y7 * z5 - y4 * z7 - y5 * z7 +
					y3 * (-z1 - z2 + z4 + z7) + y2 * (-z1 + z3 - z5 + z7));
				g[aE[j][0]][1] += 4 * (x4 * z1 + x5 * z1 + x1 * z2 + x5 * z2 - x7 * z2 + x1 * z3 - x4 * z3 - x7 * z3 -
					x1 * z4 - x5 * z4 + x7 * z4 - x1 * z5 + x4 * z5 + x7 * z5 - x4 * z7 - x5 * z7 +
					x3 * (-z1 - z2 + z4 + z7) + x2 * (-z1 + z3 - z5 + z7));
				g[aE[j][0]][2] += -4 * (x4 * y1 + x5 * y1 + x1 * y2 + x5 * y2 - x7 * y2 + x1 * y3 - x4 * y3 - x7 * y3 -
					x1 * y4 - x5 * y4 + x7 * y4 - x1 * y5 + x4 * y5 + x7 * y5 - x4 * y7 - x5 * y7 +
					x3 * (-y1 - y2 + y4 + y7) + x2 * (-y1 + y3 - y5 + y7));
				g[aE[j][1]][0] += 4 * (y4 * z0 + y5 * z0 + y0 * z2 - y5 * z2 - y6 * z2 + y0 * z3 + y4 * z3 - y6 * z3 -
					y0 * z4 + y5 * z4 + y6 * z4 - y0 * z5 - y4 * z5 + y6 * z5 - y4 * z6 - y5 * z6 +
					y3 * (-z0 + z2 - z4 + z6) + y2 * (-z0 - z3 + z5 + z6));
				g[aE[j][1]][1] += -4 * (x4 * z0 + x5 * z0 + x0 * z2 - x5 * z2 - x6 * z2 + x0 * z3 + x4 * z3 - x6 * z3 -
					x0 * z4 + x5 * z4 + x6 * z4 - x0 * z5 - x4 * z5 + x6 * z5 - x4 * z6 - x5 * z6 +
					x3 * (-z0 + z2 - z4 + z6) + x2 * (-z0 - z3 + z5 + z6));
				g[aE[j][1]][2] += 4 * (x4 * y0 + x5 * y0 + x0 * y2 - x5 * y2 - x6 * y2 + x0 * y3 + x4 * y3 - x6 * y3 -
					x0 * y4 + x5 * y4 + x6 * y4 - x0 * y5 - x4 * y5 + x6 * y5 - x4 * y6 - x5 * y6 +
					x3 * (-y0 + y2 - y4 + y6) + x2 * (-y0 - y3 + y5 + y6));
				g[aE[j][2]][0] += -4 * (-y5 * z0 + y7 * z0 + y0 * z1 - y5 * z1 - y6 * z1 - y0 * z3 + y6 * z3 + y7 * z3 +
					y0 * z5 - y6 * z5 - y7 * z5 + y5 * z6 - y7 * z6 + y1 * (-z0 - z3 + z5 + z6) +
					y3 * (z0 + z1 - z6 - z7) - y0 * z7 + y5 * z7 + y6 * z7);
				g[aE[j][2]][1] += 4 * (-x5 * z0 + x7 * z0 + x0 * z1 - x5 * z1 - x6 * z1 - x0 * z3 + x6 * z3 + x7 * z3 +
					x0 * z5 - x6 * z5 - x7 * z5 + x5 * z6 - x7 * z6 + x1 * (-z0 - z3 + z5 + z6) +
					x3 * (z0 + z1 - z6 - z7) - x0 * z7 + x5 * z7 + x6 * z7);
				g[aE[j][2]][2] += -4 * (-x5 * y0 + x7 * y0 + x0 * y1 - x5 * y1 - x6 * y1 - x0 * y3 + x6 * y3 + x7 * y3 +
					x0 * y5 - x6 * y5 - x7 * y5 + x5 * y6 - x7 * y6 + x1 * (-y0 - y3 + y5 + y6) +
					x3 * (y0 + y1 - y6 - y7) - x0 * y7 + x5 * y7 + x6 * y7);
				g[aE[j][3]][0] += -4 * (y4 * z0 + y7 * z0 + y0 * z1 + y4 * z1 - y6 * z1 + y0 * z2 - y6 * z2 - y7 * z2 -
					y0 * z4 + y6 * z4 + y7 * z4 - y4 * z6 - y7 * z6 + y1 * (-z0 + z2 - z4 + z6) -
					y0 * z7 - y4 * z7 + y6 * z7 + y2 * (-z0 - z1 + z6 + z7));
				g[aE[j][3]][1] += 4 * (x4 * z0 + x7 * z0 + x0 * z1 + x4 * z1 - x6 * z1 + x0 * z2 - x6 * z2 - x7 * z2 -
					x0 * z4 + x6 * z4 + x7 * z4 - x4 * z6 - x7 * z6 + x1 * (-z0 + z2 - z4 + z6) -
					x0 * z7 - x4 * z7 + x6 * z7 + x2 * (-z0 - z1 + z6 + z7));
				g[aE[j][3]][2] += -4 * (x4 * y0 + x7 * y0 + x0 * y1 + x4 * y1 - x6 * y1 + x0 * y2 - x6 * y2 - x7 * y2 -
					x0 * y4 + x6 * y4 + x7 * y4 - x4 * y6 - x7 * y6 + x1 * (-y0 + y2 - y4 + y6) -
					x0 * y7 - x4 * y7 + x6 * y7 + x2 * (-y0 - y1 + y6 + y7));
				g[aE[j][4]][0] += 4 * (-y5 * z0 + y7 * z0 + y0 * z1 - y5 * z1 - y6 * z1 - y0 * z3 + y6 * z3 + y7 * z3 +
					y0 * z5 - y6 * z5 - y7 * z5 + y5 * z6 - y7 * z6 + y1 * (-z0 - z3 + z5 + z6) +
					y3 * (z0 + z1 - z6 - z7) - y0 * z7 + y5 * z7 + y6 * z7);
				g[aE[j][4]][1] += -4 * (-x5 * z0 + x7 * z0 + x0 * z1 - x5 * z1 - x6 * z1 - x0 * z3 + x6 * z3 + x7 * z3 +
					x0 * z5 - x6 * z5 - x7 * z5 + x5 * z6 - x7 * z6 + x1 * (-z0 - z3 + z5 + z6) +
					x3 * (z0 + z1 - z6 - z7) - x0 * z7 + x5 * z7 + x6 * z7);
				g[aE[j][4]][2] += 4 * (-x5 * y0 + x7 * y0 + x0 * y1 - x5 * y1 - x6 * y1 - x0 * y3 + x6 * y3 + x7 * y3 +
					x0 * y5 - x6 * y5 - x7 * y5 + x5 * y6 - x7 * y6 + x1 * (-y0 - y3 + y5 + y6) +
					x3 * (y0 + y1 - y6 - y7) - x0 * y7 + x5 * y7 + x6 * y7);
				g[aE[j][5]][0] += 4 * (y4 * z0 + y7 * z0 + y0 * z1 + y4 * z1 - y6 * z1 + y0 * z2 - y6 * z2 - y7 * z2 -
					y0 * z4 + y6 * z4 + y7 * z4 - y4 * z6 - y7 * z6 + y1 * (-z0 + z2 - z4 + z6) -
					y0 * z7 - y4 * z7 + y6 * z7 + y2 * (-z0 - z1 + z6 + z7));
				g[aE[j][5]][1] += -4 * (x4 * z0 + x7 * z0 + x0 * z1 + x4 * z1 - x6 * z1 + x0 * z2 - x6 * z2 - x7 * z2 -
					x0 * z4 + x6 * z4 + x7 * z4 - x4 * z6 - x7 * z6 + x1 * (-z0 + z2 - z4 + z6) -
					x0 * z7 - x4 * z7 + x6 * z7 + x2 * (-z0 - z1 + z6 + z7));
				g[aE[j][5]][2] += 4 * (x4 * y0 + x7 * y0 + x0 * y1 + x4 * y1 - x6 * y1 + x0 * y2 - x6 * y2 - x7 * y2 -
					x0 * y4 + x6 * y4 + x7 * y4 - x4 * y6 - x7 * y6 + x1 * (-y0 + y2 - y4 + y6) -
					x0 * y7 - x4 * y7 + x6 * y7 + x2 * (-y0 - y1 + y6 + y7));
				g[aE[j][6]][0] += 4 * (y4 * z1 + y5 * z1 + y1 * z2 + y5 * z2 - y7 * z2 + y1 * z3 - y4 * z3 - y7 * z3 -
					y1 * z4 - y5 * z4 + y7 * z4 - y1 * z5 + y4 * z5 + y7 * z5 - y4 * z7 - y5 * z7 +
					y3 * (-z1 - z2 + z4 + z7) + y2 * (-z1 + z3 - z5 + z7));
				g[aE[j][6]][1] += -4 * (x4 * z1 + x5 * z1 + x1 * z2 + x5 * z2 - x7 * z2 + x1 * z3 - x4 * z3 - x7 * z3 -
					x1 * z4 - x5 * z4 + x7 * z4 - x1 * z5 + x4 * z5 + x7 * z5 - x4 * z7 - x5 * z7 +
					x3 * (-z1 - z2 + z4 + z7) + x2 * (-z1 + z3 - z5 + z7));
				g[aE[j][6]][2] += 4 * (x4 * y1 + x5 * y1 + x1 * y2 + x5 * y2 - x7 * y2 + x1 * y3 - x4 * y3 - x7 * y3 -
					x1 * y4 - x5 * y4 + x7 * y4 - x1 * y5 + x4 * y5 + x7 * y5 - x4 * y7 - x5 * y7 +
					x3 * (-y1 - y2 + y4 + y7) + x2 * (-y1 + y3 - y5 + y7));
				g[aE[j][7]][0] += -4 * (y4 * z0 + y5 * z0 + y0 * z2 - y5 * z2 - y6 * z2 + y0 * z3 + y4 * z3 - y6 * z3 -
					y0 * z4 + y5 * z4 + y6 * z4 - y0 * z5 - y4 * z5 + y6 * z5 - y4 * z6 - y5 * z6 +
					y3 * (-z0 + z2 - z4 + z6) + y2 * (-z0 - z3 + z5 + z6));
				g[aE[j][7]][1] += 4 * (x4 * z0 + x5 * z0 + x0 * z2 - x5 * z2 - x6 * z2 + x0 * z3 + x4 * z3 - x6 * z3 -
					x0 * z4 + x5 * z4 + x6 * z4 - x0 * z5 - x4 * z5 + x6 * z5 - x4 * z6 - x5 * z6 +
					x3 * (-z0 + z2 - z4 + z6) + x2 * (-z0 - z3 + z5 + z6));
				g[aE[j][7]][2] += -4 * (x4 * y0 + x5 * y0 + x0 * y2 - x5 * y2 - x6 * y2 + x0 * y3 + x4 * y3 - x6 * y3 -
					x0 * y4 + x5 * y4 + x6 * y4 - x0 * y5 - x4 * y5 + x6 * y5 - x4 * y6 - x5 * y6 +
					x3 * (-y0 + y2 - y4 + y6) + x2 * (-y0 - y3 + y5 + y6));
			}
			else if (minIdx[0] == 1) {
				pair = true;
				g[aE[j][0]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (-2 * z2 - 2 * z3 + 2 * z4 + 2 * z5) + (2 * y2 + 2 * y3 - 2 * y4 - 2 * y5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) -
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) - (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][0]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * z2 + 2 * z3 - 2 * z4 - 2 * z5) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(-2 * x2 - 2 * x3 + 2 * x4 + 2 * x5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) - (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][0]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * y2 - 2 * y3 + 2 * y4 + 2 * y5) + (2 * x2 + 2 * x3 - 2 * x4 - 2 * x5) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) -
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][1]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (-2 * z2 - 2 * z3 + 2 * z4 + 2 * z5) + (2 * y2 + 2 * y3 - 2 * y4 - 2 * y5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) +
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) + (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][1]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * z2 + 2 * z3 - 2 * z4 - 2 * z5) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(-2 * x2 - 2 * x3 + 2 * x4 + 2 * x5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][1]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * y2 - 2 * y3 + 2 * y4 + 2 * y5) + (2 * x2 + 2 * x3 - 2 * x4 - 2 * x5) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) +
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][2]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (2 * z0 + 2 * z1 - 2 * z6 - 2 * z7) + (-2 * y0 - 2 * y1 + 2 * y6 + 2 * y7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) +
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) + (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][2]][1] += ((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (2 * x0 + 2 * x1 - 2 * x6 - 2 * x7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) +
					(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) + (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * z0 - 2 * z1 + 2 * z6 + 2 * z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][2]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * y0 + 2 * y1 - 2 * y6 - 2 * y7) + (-2 * x0 - 2 * x1 + 2 * x6 + 2 * x7) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) +
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][3]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (2 * z0 + 2 * z1 - 2 * z6 - 2 * z7) + (-2 * y0 - 2 * y1 + 2 * y6 + 2 * y7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) -
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) - (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][3]][1] += (-((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7)) + (2 * x0 + 2 * x1 - 2 * x6 - 2 * x7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) -
					(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) + (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * z0 - 2 * z1 + 2 * z6 + 2 * z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][3]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * y0 + 2 * y1 - 2 * y6 - 2 * y7) + (-2 * x0 - 2 * x1 + 2 * x6 + 2 * x7) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) -
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(-2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][4]][0] += ((2 * y0 + 2 * y1 - 2 * y6 - 2 * y7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) - (y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) -
					(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (-2 * z0 - 2 * z1 + 2 * z6 + 2 * z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][4]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * z0 + 2 * z1 - 2 * z6 - 2 * z7) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(-2 * x0 - 2 * x1 + 2 * x6 + 2 * x7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) - (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][4]][2] += ((2 * x0 + 2 * x1 - 2 * x6 - 2 * x7) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) - (x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) -
					(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) + (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * y0 - 2 * y1 + 2 * y6 + 2 * y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][5]][0] += ((2 * y0 + 2 * y1 - 2 * y6 - 2 * y7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
					(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (-2 * z0 - 2 * z1 + 2 * z6 + 2 * z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][5]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * z0 + 2 * z1 - 2 * z6 - 2 * z7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(-2 * x0 - 2 * x1 + 2 * x6 + 2 * x7) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][5]][2] += ((2 * x0 + 2 * x1 - 2 * x6 - 2 * x7) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) + (x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) +
					(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) + (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * y0 - 2 * y1 + 2 * y6 + 2 * y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) - 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][6]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (2 * z2 + 2 * z3 - 2 * z4 - 2 * z5) + (-2 * y2 - 2 * y3 + 2 * y4 + 2 * y5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) +
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) + (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][6]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * z2 - 2 * z3 + 2 * z4 + 2 * z5) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(2 * x2 + 2 * x3 - 2 * x4 - 2 * x5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][6]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * y2 + 2 * y3 - 2 * y4 - 2 * y5) + (-2 * x2 - 2 * x3 + 2 * x4 + 2 * x5) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) +
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) +
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][7]][0] += ((-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (2 * z2 + 2 * z3 - 2 * z4 - 2 * z5) + (-2 * y2 - 2 * y3 + 2 * y4 + 2 * y5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) -
					(y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) - (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][7]][1] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (-2 * z2 - 2 * z3 + 2 * z4 + 2 * z5) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) +
					(2 * x2 + 2 * x3 - 2 * x4 - 2 * x5) * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) - (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7) * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
							pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) *
							(pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
				g[aE[j][7]][2] += ((-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * (2 * y2 + 2 * y3 - 2 * y4 - 2 * y5) + (-2 * x2 - 2 * x3 + 2 * x4 + 2 * x5) * (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) -
					(x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) - (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 0.5) -
					(0.5 * (((x0 + x1 + x2 + x3 - x4 - x5 - x6 - x7) * (-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) + (-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7)) *
						(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) + (-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7) *
						((-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7) * (z0 + z1 + z2 + z3 - z4 - z5 - z6 - z7) + (-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7)) +
						(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7) * ((y0 + y1 + y2 + y3 - y4 - y5 - y6 - y7) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) +
							(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7) * (-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7))) *
						(2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
							(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) + 2 * (pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) +
									pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) * (-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7) *
							(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)) -
							2 * (-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7) * (pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) +
								pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) * (pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) +
									pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)))) /
					pow((pow(-x0 + x1 + x2 - x3 - x4 + x5 + x6 - x7, 2) + pow(-y0 + y1 + y2 - y3 - y4 + y5 + y6 - y7, 2) + pow(-z0 + z1 + z2 - z3 - z4 + z5 + z6 - z7, 2)) *
						(pow(-x0 - x1 + x2 + x3 - x4 - x5 + x6 + x7, 2) + pow(-y0 - y1 + y2 + y3 - y4 - y5 + y6 + y7, 2) + pow(-z0 - z1 + z2 + z3 - z4 - z5 + z6 + z7, 2)) *
						(pow(-x0 - x1 - x2 - x3 + x4 + x5 + x6 + x7, 2) + pow(-y0 - y1 - y2 - y3 + y4 + y5 + y6 + y7, 2) + pow(-z0 - z1 - z2 - z3 + z4 + z5 + z6 + z7, 2)), 1.5);
			}
			// 1
			if (minIdx[1] == 0) {
				pair = true;
				g[aE[j][0]][0] += y4 * (-z1 + z3) + y3 * (z1 - z4) + y1 * (-z3 + z4);
				g[aE[j][0]][1] += x4 * (z1 - z3) + x1 * (z3 - z4) + x3 * (-z1 + z4);
				g[aE[j][0]][2] += x4 * (-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4);
				g[aE[j][1]][0] += y4 * (z0 - z3) + y0 * (z3 - z4) + y3 * (-z0 + z4);
				g[aE[j][1]][1] += x4 * (-z0 + z3) + x3 * (z0 - z4) + x0 * (-z3 + z4);
				g[aE[j][1]][2] += x4 * (y0 - y3) + x0 * (y3 - y4) + x3 * (-y0 + y4);
				g[aE[j][3]][0] += y4 * (-z0 + z1) + y1 * (z0 - z4) + y0 * (-z1 + z4);
				g[aE[j][3]][1] += x4 * (z0 - z1) + x0 * (z1 - z4) + x1 * (-z0 + z4);
				g[aE[j][3]][2] += x4 * (-y0 + y1) + x1 * (y0 - y4) + x0 * (-y1 + y4);
				g[aE[j][4]][0] += y3 * (z0 - z1) + y0 * (z1 - z3) + y1 * (-z0 + z3);
				g[aE[j][4]][1] += x3 * (-z0 + z1) + x1 * (z0 - z3) + x0 * (-z1 + z3);
				g[aE[j][4]][2] += x3 * (y0 - y1) + x0 * (y1 - y3) + x1 * (-y0 + y3);
			}
			else if (minIdx[1] == 1) {
				pair = true;
				g[aE[j][0]][0] += (y4 * (-z1 + z3) + y3 * (z1 - z4) + y1 * (-z3 + z4)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						0.5) - (0.5 * (-2 * (-x0 + x4) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) -
							2 * (-x0 + x3) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) -
							2 * (-x0 + x1) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) *
							((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) + (x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						1.5);
				g[aE[j][0]][1] += (x4 * (z1 - z3) + x1 * (z3 - z4) + x3 * (-z1 + z4)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						0.5) - (0.5 * (-2 * (-y0 + y4) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) -
							2 * (-y0 + y3) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) -
							2 * (-y0 + y1) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) *
							((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) + (x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						1.5);
				g[aE[j][0]][2] += (x4 * (-y1 + y3) + x3 * (y1 - y4) + x1 * (-y3 + y4)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						0.5) - (0.5 * (-2 * (-z0 + z1) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) -
							2 * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (-z0 + z3) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) -
							2 * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (-z0 + z4)) *
							((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) + (x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)),
						1.5);
				g[aE[j][1]][0] += (y4 * (z0 - z3) + y0 * (z3 - z4) + y3 * (-z0 + z4) - (1. * (-x0 + x1) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][1]][1] += (x4 * (-z0 + z3) + x3 * (z0 - z4) + x0 * (-z3 + z4) - (1. * (-y0 + y1) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][1]][2] += (x4 * (y0 - y3) + x0 * (y3 - y4) + x3 * (-y0 + y4) - (1. * (-z0 + z1) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][3]][0] += (y4 * (-z0 + z1) + y1 * (z0 - z4) + y0 * (-z1 + z4) - (1. * (-x0 + x3) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][3]][1] += (x4 * (z0 - z1) + x0 * (z1 - z4) + x1 * (-z0 + z4) - (1. * (-y0 + y3) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][3]][2] += (x4 * (-y0 + y1) + x1 * (y0 - y4) + x0 * (-y1 + y4) - (1. * (-z0 + z3) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][4]][0] += (y3 * (z0 - z1) + y0 * (z1 - z3) + y1 * (-z0 + z3) - (1. * (-x0 + x4) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][4]][1] += (x3 * (-z0 + z1) + x1 * (z0 - z3) + x0 * (-z1 + z3) - (1. * (-y0 + y4) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
				g[aE[j][4]][2] += (x3 * (y0 - y1) + x0 * (y1 - y3) + x1 * (-y0 + y3) - (1. * (-z0 + z4) * ((x4 * (-y0 + y3) + x3 * (y0 - y4) + x0 * (-y3 + y4)) * (z0 - z1) + (y0 - y1) * (x4 * (z0 - z3) + x0 * (z3 - z4) + x3 * (-z0 + z4)) +
					(x0 - x1) * (y4 * (-z0 + z3) + y3 * (z0 - z4) + y0 * (-z3 + z4)))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)), 0.5);
			}
			// 2
			if (minIdx[2] == 0) {
				pair = true;
				g[aE[j][0]][0] += y5 * (-z1 + z2) + y2 * (z1 - z5) + y1 * (-z2 + z5);
				g[aE[j][0]][1] += x5 * (z1 - z2) + x1 * (z2 - z5) + x2 * (-z1 + z5);
				g[aE[j][0]][2] += x5 * (-y1 + y2) + x2 * (y1 - y5) + x1 * (-y2 + y5);
				g[aE[j][1]][0] += y5 * (z0 - z2) + y0 * (z2 - z5) + y2 * (-z0 + z5);
				g[aE[j][1]][1] += x5 * (-z0 + z2) + x2 * (z0 - z5) + x0 * (-z2 + z5);
				g[aE[j][1]][2] += x5 * (y0 - y2) + x0 * (y2 - y5) + x2 * (-y0 + y5);
				g[aE[j][2]][0] += y5 * (-z0 + z1) + y1 * (z0 - z5) + y0 * (-z1 + z5);
				g[aE[j][2]][1] += x5 * (z0 - z1) + x0 * (z1 - z5) + x1 * (-z0 + z5);
				g[aE[j][2]][2] += x5 * (-y0 + y1) + x1 * (y0 - y5) + x0 * (-y1 + y5);
				g[aE[j][5]][0] += y2 * (z0 - z1) + y0 * (z1 - z2) + y1 * (-z0 + z2);
				g[aE[j][5]][1] += x2 * (-z0 + z1) + x1 * (z0 - z2) + x0 * (-z1 + z2);
				g[aE[j][5]][2] += x2 * (y0 - y1) + x0 * (y1 - y2) + x1 * (-y0 + y2);
			}
			else if (minIdx[2] == 1) {
				pair = true;
				g[aE[j][0]][0] += (y5 * (-z1 + z2) + y2 * (z1 - z5) + y1 * (-z2 + z5) - (1. * (x0 - x1) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][0]][1] += (x5 * (z1 - z2) + x1 * (z2 - z5) + x2 * (-z1 + z5) - (1. * (y0 - y1) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][0]][2] += (x5 * (-y1 + y2) + x2 * (y1 - y5) + x1 * (-y2 + y5) - (1. * (z0 - z1) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][1]][0] += (y5 * (z0 - z2) + y0 * (z2 - z5) + y2 * (-z0 + z5)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						0.5) - (0.5 * (-2 * (-x1 + x5) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) -
							2 * (-x1 + x2) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) -
							2 * (x0 - x1) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) *
							((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) + (y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						1.5);
				g[aE[j][1]][1] += (x5 * (-z0 + z2) + x2 * (z0 - z5) + x0 * (-z2 + z5)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						0.5) - (0.5 * (-2 * (-y1 + y5) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) -
							2 * (-y1 + y2) * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) -
							2 * (y0 - y1) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) *
							((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) + (y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						1.5);
				g[aE[j][1]][2] += (x5 * (y0 - y2) + x0 * (y2 - y5) + x2 * (-y0 + y5)) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						0.5) - (0.5 * (-2 * (z0 - z1) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) -
							2 * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (-z1 + z2) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) -
							2 * (pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (-z1 + z5)) *
							((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) + (y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)),
						1.5);
				g[aE[j][2]][0] += (y5 * (-z0 + z1) + y1 * (z0 - z5) + y0 * (-z1 + z5) - (1. * (-x1 + x2) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][2]][1] += (x5 * (z0 - z1) + x0 * (z1 - z5) + x1 * (-z0 + z5) - (1. * (-y1 + y2) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][2]][2] += (x5 * (-y0 + y1) + x1 * (y0 - y5) + x0 * (-y1 + y5) - (1. * (-z1 + z2) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][5]][0] += (y2 * (z0 - z1) + y0 * (z1 - z2) + y1 * (-z0 + z2) - (1. * (-x1 + x5) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][5]][1] += (x2 * (-z0 + z1) + x1 * (z0 - z2) + x0 * (-z1 + z2) - (1. * (-y1 + y5) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
				g[aE[j][5]][2] += (x2 * (y0 - y1) + x0 * (y1 - y2) + x1 * (-y0 + y2) - (1. * (-z1 + z5) * ((x5 * (y0 - y1) + x0 * (y1 - y5) + x1 * (-y0 + y5)) * (z1 - z2) + (x1 - x2) * (y5 * (z0 - z1) + y0 * (z1 - z5) + y1 * (-z0 + z5)) +
					(y1 - y2) * (x5 * (-z0 + z1) + x1 * (z0 - z5) + x0 * (-z1 + z5)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x0 - x1, 2) + pow(y0 - y1, 2) + pow(z0 - z1, 2)) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)), 0.5);
			}
			// 3
			if (minIdx[3] == 0) {
				pair = true;
				g[aE[j][1]][0] += y6 * (-z2 + z3) + y3 * (z2 - z6) + y2 * (-z3 + z6);
				g[aE[j][1]][1] += x6 * (z2 - z3) + x2 * (z3 - z6) + x3 * (-z2 + z6);
				g[aE[j][1]][2] += x6 * (-y2 + y3) + x3 * (y2 - y6) + x2 * (-y3 + y6);
				g[aE[j][2]][0] += y6 * (z1 - z3) + y1 * (z3 - z6) + y3 * (-z1 + z6);
				g[aE[j][2]][1] += x6 * (-z1 + z3) + x3 * (z1 - z6) + x1 * (-z3 + z6);
				g[aE[j][2]][2] += x6 * (y1 - y3) + x1 * (y3 - y6) + x3 * (-y1 + y6);
				g[aE[j][3]][0] += y6 * (-z1 + z2) + y2 * (z1 - z6) + y1 * (-z2 + z6);
				g[aE[j][3]][1] += x6 * (z1 - z2) + x1 * (z2 - z6) + x2 * (-z1 + z6);
				g[aE[j][3]][2] += x6 * (-y1 + y2) + x2 * (y1 - y6) + x1 * (-y2 + y6);
				g[aE[j][6]][0] += y3 * (z1 - z2) + y1 * (z2 - z3) + y2 * (-z1 + z3);
				g[aE[j][6]][1] += x3 * (-z1 + z2) + x2 * (z1 - z3) + x1 * (-z2 + z3);
				g[aE[j][6]][2] += x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3);
			}
			else if (minIdx[3] == 1) {
				pair = true;
				g[aE[j][1]][0] += (y6 * (-z2 + z3) + y3 * (z2 - z6) + y2 * (-z3 + z6) - (1. * (x1 - x2) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][1]][1] += (x6 * (z2 - z3) + x2 * (z3 - z6) + x3 * (-z2 + z6) - (1. * (y1 - y2) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][1]][2] += (x6 * (-y2 + y3) + x3 * (y2 - y6) + x2 * (-y3 + y6) - (1. * (z1 - z2) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][2]][0] += (y6 * (z1 - z3) + y1 * (z3 - z6) + y3 * (-z1 + z6)) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						0.5) - (0.5 * (-2 * (-x2 + x6) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) -
							2 * (-x2 + x3) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) -
							2 * (x1 - x2) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) *
							((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) + (y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						1.5);
				g[aE[j][2]][1] += (x6 * (-z1 + z3) + x3 * (z1 - z6) + x1 * (-z3 + z6)) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						0.5) - (0.5 * (-2 * (-y2 + y6) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) -
							2 * (-y2 + y3) * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) -
							2 * (y1 - y2) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) *
							((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) + (y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						1.5);
				g[aE[j][2]][2] += (x6 * (y1 - y3) + x1 * (y3 - y6) + x3 * (-y1 + y6)) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						0.5) - (0.5 * (-2 * (z1 - z2) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) -
							2 * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (-z2 + z3) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) -
							2 * (pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (-z2 + z6)) *
							((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) + (y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)),
						1.5);
				g[aE[j][3]][0] += (y6 * (-z1 + z2) + y2 * (z1 - z6) + y1 * (-z2 + z6) - (1. * (-x2 + x3) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][3]][1] += (x6 * (z1 - z2) + x1 * (z2 - z6) + x2 * (-z1 + z6) - (1. * (-y2 + y3) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][3]][2] += (x6 * (-y1 + y2) + x2 * (y1 - y6) + x1 * (-y2 + y6) - (1. * (-z2 + z3) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][6]][0] += (y3 * (z1 - z2) + y1 * (z2 - z3) + y2 * (-z1 + z3) - (1. * (-x2 + x6) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][6]][1] += (x3 * (-z1 + z2) + x2 * (z1 - z3) + x1 * (-z2 + z3) - (1. * (-y2 + y6) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
				g[aE[j][6]][2] += (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3) - (1. * (-z2 + z6) * ((x6 * (y1 - y2) + x1 * (y2 - y6) + x2 * (-y1 + y6)) * (z2 - z3) + (x2 - x3) * (y6 * (z1 - z2) + y1 * (z2 - z6) + y2 * (-z1 + z6)) +
					(y2 - y3) * (x6 * (-z1 + z2) + x2 * (z1 - z6) + x1 * (-z2 + z6)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)), 0.5);
			}
			// 4
			if (minIdx[4] == 0) {
				pair = true;
				g[aE[j][0]][0] += y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7);
				g[aE[j][0]][1] += x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7);
				g[aE[j][0]][2] += x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7);
				g[aE[j][2]][0] += y7 * (z0 - z3) + y0 * (z3 - z7) + y3 * (-z0 + z7);
				g[aE[j][2]][1] += x7 * (-z0 + z3) + x3 * (z0 - z7) + x0 * (-z3 + z7);
				g[aE[j][2]][2] += x7 * (y0 - y3) + x0 * (y3 - y7) + x3 * (-y0 + y7);
				g[aE[j][3]][0] += y7 * (-z0 + z2) + y2 * (z0 - z7) + y0 * (-z2 + z7);
				g[aE[j][3]][1] += x7 * (z0 - z2) + x0 * (z2 - z7) + x2 * (-z0 + z7);
				g[aE[j][3]][2] += x7 * (-y0 + y2) + x2 * (y0 - y7) + x0 * (-y2 + y7);
				g[aE[j][7]][0] += y3 * (z0 - z2) + y0 * (z2 - z3) + y2 * (-z0 + z3);
				g[aE[j][7]][1] += x3 * (-z0 + z2) + x2 * (z0 - z3) + x0 * (-z2 + z3);
				g[aE[j][7]][2] += x3 * (y0 - y2) + x0 * (y2 - y3) + x2 * (-y0 + y3);
			}
			else if (minIdx[4] == 1) {
				pair = true;
				g[aE[j][0]][0] += (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7) - (1. * (x0 - x3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][0]][1] += (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7) - (1. * (y0 - y3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][0]][2] += (x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7) - (1. * (z0 - z3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][2]][0] += (y7 * (z0 - z3) + y0 * (z3 - z7) + y3 * (-z0 + z7) - (1. * (x2 - x3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][2]][1] += (x7 * (-z0 + z3) + x3 * (z0 - z7) + x0 * (-z3 + z7) - (1. * (y2 - y3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][2]][2] += (x7 * (y0 - y3) + x0 * (y3 - y7) + x3 * (-y0 + y7) - (1. * (z2 - z3) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][3]][0] += (y7 * (-z0 + z2) + y2 * (z0 - z7) + y0 * (-z2 + z7)) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						0.5) - (0.5 * (-2 * (-x3 + x7) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) -
							2 * (x2 - x3) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) -
							2 * (x0 - x3) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) *
							((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) + (x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						1.5);
				g[aE[j][3]][1] += (x7 * (z0 - z2) + x0 * (z2 - z7) + x2 * (-z0 + z7)) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						0.5) - (0.5 * (-2 * (-y3 + y7) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) -
							2 * (y2 - y3) * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) -
							2 * (y0 - y3) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) *
							((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) + (x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						1.5);
				g[aE[j][3]][2] += (x7 * (-y0 + y2) + x2 * (y0 - y7) + x0 * (-y2 + y7)) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						0.5) - (0.5 * (-2 * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (z0 - z3) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) -
							2 * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (z2 - z3) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) -
							2 * (pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (-z3 + z7)) *
							((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) + (x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)),
						1.5);
				g[aE[j][7]][0] += (y3 * (z0 - z2) + y0 * (z2 - z3) + y2 * (-z0 + z3) - (1. * (-x3 + x7) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][7]][1] += (x3 * (-z0 + z2) + x2 * (z0 - z3) + x0 * (-z2 + z3) - (1. * (-y3 + y7) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
				g[aE[j][7]][2] += (x3 * (y0 - y2) + x0 * (y2 - y3) + x2 * (-y0 + y3) - (1. * (-z3 + z7) * ((x7 * (-y2 + y3) + x3 * (y2 - y7) + x2 * (-y3 + y7)) * (z0 - z3) + (y0 - y3) * (x7 * (z2 - z3) + x2 * (z3 - z7) + x3 * (-z2 + z7)) +
					(x0 - x3) * (y7 * (-z2 + z3) + y3 * (z2 - z7) + y2 * (-z3 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x0 - x3, 2) + pow(y0 - y3, 2) + pow(z0 - z3, 2)) * (pow(x2 - x3, 2) + pow(y2 - y3, 2) + pow(z2 - z3, 2)) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)), 0.5);
			}
			// 5
			if (minIdx[5] == 0) {
				pair = true;
				g[aE[j][0]][0] += y7 * (-z4 + z5) + y5 * (z4 - z7) + y4 * (-z5 + z7);
				g[aE[j][0]][1] += x7 * (z4 - z5) + x4 * (z5 - z7) + x5 * (-z4 + z7);
				g[aE[j][0]][2] += x7 * (-y4 + y5) + x5 * (y4 - y7) + x4 * (-y5 + y7);
				g[aE[j][4]][0] += y7 * (z0 - z5) + y0 * (z5 - z7) + y5 * (-z0 + z7);
				g[aE[j][4]][1] += x7 * (-z0 + z5) + x5 * (z0 - z7) + x0 * (-z5 + z7);
				g[aE[j][4]][2] += x7 * (y0 - y5) + x0 * (y5 - y7) + x5 * (-y0 + y7);
				g[aE[j][5]][0] += y7 * (-z0 + z4) + y4 * (z0 - z7) + y0 * (-z4 + z7);
				g[aE[j][5]][1] += x7 * (z0 - z4) + x0 * (z4 - z7) + x4 * (-z0 + z7);
				g[aE[j][5]][2] += x7 * (-y0 + y4) + x4 * (y0 - y7) + x0 * (-y4 + y7);
				g[aE[j][7]][0] += y5 * (z0 - z4) + y0 * (z4 - z5) + y4 * (-z0 + z5);
				g[aE[j][7]][1] += x5 * (-z0 + z4) + x4 * (z0 - z5) + x0 * (-z4 + z5);
				g[aE[j][7]][2] += x5 * (y0 - y4) + x0 * (y4 - y5) + x4 * (-y0 + y5);
			}
			else if (minIdx[5] == 1) {
				pair = true;
				g[aE[j][0]][0] += (y7 * (-z4 + z5) - (1. * (x0 - x4) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) + y5 * (z4 - z7) + y4 * (-z5 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][0]][1] += (x7 * (z4 - z5) - (1. * (y0 - y4) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) + x4 * (z5 - z7) + x5 * (-z4 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][0]][2] += (x7 * (-y4 + y5) + x5 * (y4 - y7) + x4 * (-y5 + y7) - (1. * (z0 - z4) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][4]][0] += (-0.5 * (-2 * (-x4 + x7) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) -
					2 * (-x4 + x5) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) -
					2 * (x0 - x4) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) *
					((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) + (x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)),
						1.5) + (y7 * (z0 - z5) + y0 * (z5 - z7) + y5 * (-z0 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)),
						0.5);
				g[aE[j][4]][1] += (-0.5 * (-2 * (-y4 + y7) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) -
					2 * (-y4 + y5) * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) -
					2 * (y0 - y4) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) *
					((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) + (x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)),
						1.5) + (x7 * (-z0 + z5) + x5 * (z0 - z7) + x0 * (-z5 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)),
						0.5);
				g[aE[j][4]][2] += (x7 * (y0 - y5) + x0 * (y5 - y7) + x5 * (-y0 + y7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)),
						0.5) - (0.5 * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
							(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7)) * (-2 * (z0 - z4) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) *
								(pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) - 2 * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (-z4 + z5) *
								(pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) - 2 * (pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) *
								(-z4 + z7))) / pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) *
									(pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 1.5);
				g[aE[j][5]][0] += (y7 * (-z0 + z4) - (1. * (-x4 + x5) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) + y4 * (z0 - z7) + y0 * (-z4 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][5]][1] += (x7 * (z0 - z4) - (1. * (-y4 + y5) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) + x0 * (z4 - z7) + x4 * (-z0 + z7)) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][5]][2] += (x7 * (-y0 + y4) + x4 * (y0 - y7) + x0 * (-y4 + y7) - (1. * (-z4 + z5) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][7]][0] += (y5 * (z0 - z4) + y0 * (z4 - z5) + y4 * (-z0 + z5) - (1. * (-x4 + x7) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][7]][1] += (x5 * (-z0 + z4) + x4 * (z0 - z5) + x0 * (-z4 + z5) - (1. * (-y4 + y7) * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7))) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
				g[aE[j][7]][2] += (x5 * (y0 - y4) + x0 * (y4 - y5) + x4 * (-y0 + y5) - (1. * ((y4 - y7) * (x5 * (z0 - z4) + x0 * (z4 - z5) + x4 * (-z0 + z5)) + (x4 - x7) * (y5 * (-z0 + z4) + y4 * (z0 - z5) + y0 * (-z4 + z5)) +
					(x5 * (-y0 + y4) + x4 * (y0 - y5) + x0 * (-y4 + y5)) * (z4 - z7)) * (-z4 + z7)) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x0 - x4, 2) + pow(y0 - y4, 2) + pow(z0 - z4, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)), 0.5);
			}
			// 6
			if (minIdx[6] == 0) {
				pair = true;
				g[aE[j][1]][0] += y6 * (-z4 + z5) + y5 * (z4 - z6) + y4 * (-z5 + z6);
				g[aE[j][1]][1] += x6 * (z4 - z5) + x4 * (z5 - z6) + x5 * (-z4 + z6);
				g[aE[j][1]][2] += x6 * (-y4 + y5) + x5 * (y4 - y6) + x4 * (-y5 + y6);
				g[aE[j][4]][0] += y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6);
				g[aE[j][4]][1] += x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6);
				g[aE[j][4]][2] += x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6);
				g[aE[j][5]][0] += y6 * (-z1 + z4) + y4 * (z1 - z6) + y1 * (-z4 + z6);
				g[aE[j][5]][1] += x6 * (z1 - z4) + x1 * (z4 - z6) + x4 * (-z1 + z6);
				g[aE[j][5]][2] += x6 * (-y1 + y4) + x4 * (y1 - y6) + x1 * (-y4 + y6);
				g[aE[j][6]][0] += y5 * (z1 - z4) + y1 * (z4 - z5) + y4 * (-z1 + z5);
				g[aE[j][6]][1] += x5 * (-z1 + z4) + x4 * (z1 - z5) + x1 * (-z4 + z5);
				g[aE[j][6]][2] += x5 * (y1 - y4) + x1 * (y4 - y5) + x4 * (-y1 + y5);
			}
			else if (minIdx[6] == 1) {
				pair = true;
				g[aE[j][1]][0] += (y6 * (-z4 + z5) + y5 * (z4 - z6) + y4 * (-z5 + z6) - (1. * (x1 - x5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][1]][1] += (x6 * (z4 - z5) + x4 * (z5 - z6) + x5 * (-z4 + z6) - (1. * (y1 - y5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][1]][2] += (x6 * (-y4 + y5) + x5 * (y4 - y6) + x4 * (-y5 + y6) - (1. * (z1 - z5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][4]][0] += (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6) - (1. * (x4 - x5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][4]][1] += (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6) - (1. * (y4 - y5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][4]][2] += (x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6) - (1. * (z4 - z5) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][5]][0] += (y6 * (-z1 + z4) + y4 * (z1 - z6) + y1 * (-z4 + z6)) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						0.5) - (0.5 * (-2 * (-x5 + x6) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) -
							2 * (x4 - x5) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (x1 - x5) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) *
							((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) + (y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						1.5);
				g[aE[j][5]][1] += (x6 * (z1 - z4) + x1 * (z4 - z6) + x4 * (-z1 + z6)) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						0.5) - (0.5 * (-2 * (-y5 + y6) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) -
							2 * (y4 - y5) * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (y1 - y5) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) *
							((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) + (y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						1.5);
				g[aE[j][5]][2] += (x6 * (-y1 + y4) + x4 * (y1 - y6) + x1 * (-y4 + y6)) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						0.5) - (0.5 * (-2 * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (z1 - z5) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (z4 - z5) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (-z5 + z6)) *
							((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) + (y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)),
						1.5);
				g[aE[j][6]][0] += (y5 * (z1 - z4) + y1 * (z4 - z5) + y4 * (-z1 + z5) - (1. * (-x5 + x6) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][6]][1] += (x5 * (-z1 + z4) + x4 * (z1 - z5) + x1 * (-z4 + z5) - (1. * (-y5 + y6) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
				g[aE[j][6]][2] += (x5 * (y1 - y4) + x1 * (y4 - y5) + x4 * (-y1 + y5) - (1. * (-z5 + z6) * ((x6 * (y1 - y5) + x1 * (y5 - y6) + x5 * (-y1 + y6)) * (z4 - z5) + (x4 - x5) * (y6 * (z1 - z5) + y1 * (z5 - z6) + y5 * (-z1 + z6)) +
					(y4 - y5) * (x6 * (-z1 + z5) + x5 * (z1 - z6) + x1 * (-z5 + z6)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x1 - x5, 2) + pow(y1 - y5, 2) + pow(z1 - z5, 2)) * (pow(x4 - x5, 2) + pow(y4 - y5, 2) + pow(z4 - z5, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)), 0.5);
			}
			// 7
			if (minIdx[7] == 0) {
				pair = true;
				g[aE[j][2]][0] += y7 * (-z5 + z6) + y6 * (z5 - z7) + y5 * (-z6 + z7);
				g[aE[j][2]][1] += x7 * (z5 - z6) + x5 * (z6 - z7) + x6 * (-z5 + z7);
				g[aE[j][2]][2] += x7 * (-y5 + y6) + x6 * (y5 - y7) + x5 * (-y6 + y7);
				g[aE[j][5]][0] += y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7);
				g[aE[j][5]][1] += x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7);
				g[aE[j][5]][2] += x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7);
				g[aE[j][6]][0] += y7 * (-z2 + z5) + y5 * (z2 - z7) + y2 * (-z5 + z7);
				g[aE[j][6]][1] += x7 * (z2 - z5) + x2 * (z5 - z7) + x5 * (-z2 + z7);
				g[aE[j][6]][2] += x7 * (-y2 + y5) + x5 * (y2 - y7) + x2 * (-y5 + y7);
				g[aE[j][7]][0] += y6 * (z2 - z5) + y2 * (z5 - z6) + y5 * (-z2 + z6);
				g[aE[j][7]][1] += x6 * (-z2 + z5) + x5 * (z2 - z6) + x2 * (-z5 + z6);
				g[aE[j][7]][2] += x6 * (y2 - y5) + x2 * (y5 - y6) + x5 * (-y2 + y6);
			}
			else if (minIdx[7] == 1) {
				pair = true;
				g[aE[j][2]][0] += (y7 * (-z5 + z6) + y6 * (z5 - z7) + y5 * (-z6 + z7) - (1. * (x2 - x6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][2]][1] += (x7 * (z5 - z6) + x5 * (z6 - z7) + x6 * (-z5 + z7) - (1. * (y2 - y6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][2]][2] += (x7 * (-y5 + y6) + x6 * (y5 - y7) + x5 * (-y6 + y7) - (1. * (z2 - z6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][5]][0] += (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7) - (1. * (x5 - x6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][5]][1] += (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7) - (1. * (y5 - y6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][5]][2] += (x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7) - (1. * (z5 - z6) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][6]][0] += (y7 * (-z2 + z5) + y5 * (z2 - z7) + y2 * (-z5 + z7)) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (-x6 + x7) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (x5 - x6) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (x2 - x6) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) *
							((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) + (y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
				g[aE[j][6]][1] += (x7 * (z2 - z5) + x2 * (z5 - z7) + x5 * (-z2 + z7)) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (-y6 + y7) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) -
							2 * (y5 - y6) * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (y2 - y6) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) *
							((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) + (y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
				g[aE[j][6]][2] += (x7 * (-y2 + y5) + x5 * (y2 - y7) + x2 * (-y5 + y7)) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (z2 - z6) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (z5 - z6) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (-z6 + z7)) *
							((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) + (y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
				g[aE[j][7]][0] += (y6 * (z2 - z5) + y2 * (z5 - z6) + y5 * (-z2 + z6) - (1. * (-x6 + x7) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][7]][1] += (x6 * (-z2 + z5) + x5 * (z2 - z6) + x2 * (-z5 + z6) - (1. * (-y6 + y7) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][7]][2] += (x6 * (y2 - y5) + x2 * (y5 - y6) + x5 * (-y2 + y6) - (1. * (-z6 + z7) * ((x7 * (y2 - y6) + x2 * (y6 - y7) + x6 * (-y2 + y7)) * (z5 - z6) + (x5 - x6) * (y7 * (z2 - z6) + y2 * (z6 - z7) + y6 * (-z2 + z7)) +
					(y5 - y6) * (x7 * (-z2 + z6) + x6 * (z2 - z7) + x2 * (-z6 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x2 - x6, 2) + pow(y2 - y6, 2) + pow(z2 - z6, 2)) * (pow(x5 - x6, 2) + pow(y5 - y6, 2) + pow(z5 - z6, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
			}
			// 8
			if (minIdx[8] == 0) {
				pair = true;
				g[aE[j][3]][0] += y7 * (-z4 + z6) + y6 * (z4 - z7) + y4 * (-z6 + z7);
				g[aE[j][3]][1] += x7 * (z4 - z6) + x4 * (z6 - z7) + x6 * (-z4 + z7);
				g[aE[j][3]][2] += x7 * (-y4 + y6) + x6 * (y4 - y7) + x4 * (-y6 + y7);
				g[aE[j][4]][0] += y7 * (z3 - z6) + y3 * (z6 - z7) + y6 * (-z3 + z7);
				g[aE[j][4]][1] += x7 * (-z3 + z6) + x6 * (z3 - z7) + x3 * (-z6 + z7);
				g[aE[j][4]][2] += x7 * (y3 - y6) + x3 * (y6 - y7) + x6 * (-y3 + y7);
				g[aE[j][6]][0] += y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7);
				g[aE[j][6]][1] += x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7);
				g[aE[j][6]][2] += x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7);
				g[aE[j][7]][0] += y6 * (z3 - z4) + y3 * (z4 - z6) + y4 * (-z3 + z6);
				g[aE[j][7]][1] += x6 * (-z3 + z4) + x4 * (z3 - z6) + x3 * (-z4 + z6);
				g[aE[j][7]][2] += x6 * (y3 - y4) + x3 * (y4 - y6) + x4 * (-y3 + y6);
			}
			else if (minIdx[8] == 1) {
				pair = true;
				g[aE[j][3]][0] += (y7 * (-z4 + z6) + y6 * (z4 - z7) + y4 * (-z6 + z7) - (1. * (x3 - x7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][3]][1] += (x7 * (z4 - z6) + x4 * (z6 - z7) + x6 * (-z4 + z7) - (1. * (y3 - y7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][3]][2] += (x7 * (-y4 + y6) + x6 * (y4 - y7) + x4 * (-y6 + y7) - (1. * (z3 - z7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][4]][0] += (y7 * (z3 - z6) + y3 * (z6 - z7) + y6 * (-z3 + z7) - (1. * (x4 - x7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][4]][1] += (x7 * (-z3 + z6) + x6 * (z3 - z7) + x3 * (-z6 + z7) - (1. * (y4 - y7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][4]][2] += (x7 * (y3 - y6) + x3 * (y6 - y7) + x6 * (-y3 + y7) - (1. * (z4 - z7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][6]][0] += (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7) - (1. * (x6 - x7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][6]][1] += (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7) - (1. * (y6 - y7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][6]][2] += (x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7) - (1. * (z6 - z7) * ((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) +
					(x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) / (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)), 0.5);
				g[aE[j][7]][0] += (y6 * (z3 - z4) + y3 * (z4 - z6) + y4 * (-z3 + z6)) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (x6 - x7) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) -
							2 * (x4 - x7) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (x3 - x7) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) *
							((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) + (x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
				g[aE[j][7]][1] += (x6 * (-z3 + z4) + x4 * (z3 - z6) + x3 * (-z4 + z6)) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (y6 - y7) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) -
							2 * (y4 - y7) * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) -
							2 * (y3 - y7) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2))) *
							((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) + (x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
				g[aE[j][7]][2] += (x6 * (y3 - y4) + x3 * (y4 - y6) + x4 * (-y3 + y6)) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						0.5) - (0.5 * (-2 * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) * (z3 - z7) -
							2 * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)) * (z4 - z7) -
							2 * (pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (z6 - z7)) *
							((x7 * (-y3 + y4) + x4 * (y3 - y7) + x3 * (-y4 + y7)) * (z6 - z7) + (y6 - y7) * (x7 * (z3 - z4) + x3 * (z4 - z7) + x4 * (-z3 + z7)) + (x6 - x7) * (y7 * (-z3 + z4) + y4 * (z3 - z7) + y3 * (-z4 + z7)))) /
					pow((pow(x3 - x7, 2) + pow(y3 - y7, 2) + pow(z3 - z7, 2)) * (pow(x4 - x7, 2) + pow(y4 - y7, 2) + pow(z4 - z7, 2)) * (pow(x6 - x7, 2) + pow(y6 - y7, 2) + pow(z6 - z7, 2)),
						1.5);
			}
			if (!pair) sJL++;
		}
		if (sJL == affElemNum) allPositive = true;
		else allPositive = false;

		// optimize point distance to surface
		if (allPositive) {
			smallDist = -1;
			for (j = sPIdx; j < 2 * sPIdx; j++) {
				dis = PointToTri(triMesh.v[triMesh.e[triNum[j - sPIdx]][0]], triMesh.v[triMesh.e[triNum[j - sPIdx]][1]],
					triMesh.v[triMesh.e[triNum[j - sPIdx]][2]], octreeMesh.v[bP[j]], target, MAX_NUM2);
				if (dis > smallDist) {
					smallDist = dis;
					maxDistIdx = bP[j];
				}
				// calculate gradient
				g[j][0] = -3 * (octreeMesh.v[bP[j]][0] - target[0]);
				g[j][1] = -3 * (octreeMesh.v[bP[j]][1] - target[1]);
				g[j][2] = -3 * (octreeMesh.v[bP[j]][2] - target[2]);
			}
		}

		// update with gradient
		for (j = 0; j < 2 * sPIdx; j++) {
			octreeMesh.v[bP[j]][0] += LEARNING_RATE * g[j][0];
			octreeMesh.v[bP[j]][1] += LEARNING_RATE * g[j][1];
			octreeMesh.v[bP[j]][2] += LEARNING_RATE * g[j][2];
		}
		for (j = 2 * sPIdx; j < 2 * sPIdx + bP2.size(); j++) {
			octreeMesh.v[bP2[j - 2 * sPIdx]][0] += LEARNING_RATE * g[j][0];
			octreeMesh.v[bP2[j - 2 * sPIdx]][1] += LEARNING_RATE * g[j][1];
			octreeMesh.v[bP2[j - 2 * sPIdx]][2] += LEARNING_RATE * g[j][2];
		}

		if (i % UPDATE_EVERY == 0) {
			for (i = 0; i < SMOOTH_EPOCH; i++) {
				if (ELEM_THRES == 0.01) {
					for (j = 0; j < sPIdx; j++) {
						tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
						for (k = 0; k < cP2[j].size(); k++) {
							tmp[0] += octreeMesh.v[cP2[j][k]][0];
							tmp[1] += octreeMesh.v[cP2[j][k]][1];
							tmp[2] += octreeMesh.v[cP2[j][k]][2];
						}
						tmp[0] /= cP2[j].size(); tmp[1] /= cP2[j].size(); tmp[2] /= cP2[j].size();
						prop = rand() / (RAND_MAX * 1.0f);
						target[0] = prop * tmp[0] + (1 - prop) * octreeMesh.v[bP[j]][0];
						target[1] = prop * tmp[1] + (1 - prop) * octreeMesh.v[bP[j]][1];
						target[2] = prop * tmp[2] + (1 - prop) * octreeMesh.v[bP[j]][2];

						pair = true;
						for (k = 0; k < cE2[j].size(); k++) {
							for (l = 0; l < 8; l++) {
								if (bP[j] == octreeMesh.e[cE2[j][k]][l]) {
									testP[l][0] = target[0]; testP[l][1] = target[1]; testP[l][2] = target[2];
								}
								else {
									testP[l][0] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][0];
									testP[l][1] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][1];
									testP[l][2] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][2];
								}
							}
							if (Sj(testP[0], testP[1], testP[2], testP[3], testP[4], testP[5], testP[6], testP[7]) <= ELEM_THRES) {
								pair = false;
								break;
							}
						}
						if (pair) {
							octreeMesh.v[bP[j]][0] = target[0];
							octreeMesh.v[bP[j]][1] = target[1];
							octreeMesh.v[bP[j]][2] = target[2];
						}
					}

					for (j = 0; j < sPIdx; j++) {
						tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
						for (k = 0; k < cP[j].size(); k++) {
							tmp[0] += octreeMesh.v[cP[j][k]][0];
							tmp[1] += octreeMesh.v[cP[j][k]][1];
							tmp[2] += octreeMesh.v[cP[j][k]][2];
						}
						tmp[0] /= cP[j].size(); tmp[1] /= cP[j].size(); tmp[2] /= cP[j].size();
						if (bP[j + sPIdx] != maxDistIdx) {
							prop = rand() / (RAND_MAX * 1.0f);
							tmp[0] = prop * tmp[0] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][0];
							tmp[1] = prop * tmp[1] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][1];
							tmp[2] = prop * tmp[2] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][2];
						}

						minDist = MAX_NUM2;
						for (k = 0; k < triMesh.eNum; k++) {
							dis = PointToTri(triMesh.v[triMesh.e[k][0]], triMesh.v[triMesh.e[k][1]], triMesh.v[triMesh.e[k][2]],
								tmp, tmp2, minDist);
							if (dis < minDist) {
								minDist = dis;
								target[0] = tmp2[0]; target[1] = tmp2[1]; target[2] = tmp2[2];
							}
						}

						if (bP[j + sPIdx] == maxDistIdx) {
							if (dragCount[j] % 10 == 2) {
								octreeMesh.v[bP[j + sPIdx]][0] = target[0];
								octreeMesh.v[bP[j + sPIdx]][1] = target[1];
								octreeMesh.v[bP[j + sPIdx]][2] = target[2];
							}
							dragCount[j]++;
						}
						pair = true;
						for (k = 0; k < cE[j].size(); k++) {
							for (l = 0; l < 8; l++) {
								if (bP[j + sPIdx] == octreeMesh.e[cE[j][k]][l]) {
									testP[l][0] = target[0]; testP[l][1] = target[1]; testP[l][2] = target[2];
								}
								else {
									testP[l][0] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][0];
									testP[l][1] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][1];
									testP[l][2] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][2];
								}
							}
							if (Sj(testP[0], testP[1], testP[2], testP[3], testP[4], testP[5], testP[6], testP[7]) <= ELEM_THRES) {
								pair = false;
								break;
							}
						}
						if (pair) {
							octreeMesh.v[bP[j + sPIdx]][0] = target[0];
							octreeMesh.v[bP[j + sPIdx]][1] = target[1];
							octreeMesh.v[bP[j + sPIdx]][2] = target[2];
						}
					}
				}
				else {
					for (j = 0; j < sPIdx; j++)
						if (bP[j + sPIdx] == maxDistIdx) {
							minDist = MAX_NUM2;
							for (k = 0; k < triMesh.eNum; k++) {
								dis = PointToTri(triMesh.v[triMesh.e[k][0]], triMesh.v[triMesh.e[k][1]], triMesh.v[triMesh.e[k][2]],
									octreeMesh.v[bP[j + sPIdx]], tmp2, minDist);
								if (dis < minDist) {
									minDist = dis;
									target[0] = tmp2[0]; target[1] = tmp2[1]; target[2] = tmp2[2];
								}
							}

							if (dragCount[j] % 10 == 2) {
								octreeMesh.v[bP[j + sPIdx]][0] = target[0];
								octreeMesh.v[bP[j + sPIdx]][1] = target[1];
								octreeMesh.v[bP[j + sPIdx]][2] = target[2];
							}
							dragCount[j]++;
						}
				}
			}
			i = 0;
			for (j = 0; j < sPIdx; j++) {
				minDist = MAX_NUM2;
				for (k = 0; k < triMesh.eNum; k++) {
					dis = PointToTri(triMesh.v[triMesh.e[k][0]], triMesh.v[triMesh.e[k][1]], triMesh.v[triMesh.e[k][2]],
						octreeMesh.v[bP[j + sPIdx]], target, minDist);
					if (dis < minDist) {
						minDist = dis;
						triNum[j] = k;
						tmp[0] = target[0]; tmp[1] = target[1]; tmp[2] = target[2];
					}
				}
				if (smallDist < DIST_THRES2) {
					octreeMesh.v[bP[j + sPIdx]][0] = tmp[0];
					octreeMesh.v[bP[j + sPIdx]][1] = tmp[1];
					octreeMesh.v[bP[j + sPIdx]][2] = tmp[2];
				}
			}

			aveDist = 0; smallDist = -1;
			for (j = sPIdx; j < 2 * sPIdx; j++) {
				dis = PointToTri(triMesh.v[triMesh.e[triNum[j - sPIdx]][0]], triMesh.v[triMesh.e[triNum[j - sPIdx]][1]],
					triMesh.v[triMesh.e[triNum[j - sPIdx]][2]], octreeMesh.v[bP[j]], target, MAX_NUM2);
				aveDist += dis;
				if (dis > smallDist) {
					smallDist = dis; maxDistIdx = bP[j];
				}
			}
			k = 0;
			for (j = 0; j < affElemNum; j++) {
				// calculate bad jacobian number
				if (Sj(octreeMesh.v[octreeMesh.e[j][0]], octreeMesh.v[octreeMesh.e[j][1]],
					octreeMesh.v[octreeMesh.e[j][2]], octreeMesh.v[octreeMesh.e[j][3]],
					octreeMesh.v[octreeMesh.e[j][4]], octreeMesh.v[octreeMesh.e[j][5]],
					octreeMesh.v[octreeMesh.e[j][6]], octreeMesh.v[octreeMesh.e[j][7]]) <= ELEM_THRES) k++;
			}
			std::cout << "badElem: " << k << " aveDist: " << aveDist / sPIdx << " maxDist: " << smallDist << " maxDistIdx: " << maxDistIdx << std::endl;
			octreeMesh.WriteToVtk(fileName, BOX_LENGTH_RATIO, START_POINT);
			if (k == 0 && smallDist < DIST_THRES2) {
				if (ELEM_THRES == 0.01) ELEM_THRES = 0.53;
				else ELEM_THRES += 0.01;

				for (j = 0; j < sPIdx; j++) {
					tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
					for (k = 0; k < cP2[j].size(); k++) {
						tmp[0] += octreeMesh.v[cP2[j][k]][0];
						tmp[1] += octreeMesh.v[cP2[j][k]][1];
						tmp[2] += octreeMesh.v[cP2[j][k]][2];
					}
					tmp[0] /= cP2[j].size(); tmp[1] /= cP2[j].size(); tmp[2] /= cP2[j].size();
					prop = rand() / (RAND_MAX * 1.0f);
					target[0] = prop * tmp[0] + (1 - prop) * octreeMesh.v[bP[j]][0];
					target[1] = prop * tmp[1] + (1 - prop) * octreeMesh.v[bP[j]][1];
					target[2] = prop * tmp[2] + (1 - prop) * octreeMesh.v[bP[j]][2];

					pair = true;
					for (k = 0; k < cE2[j].size(); k++) {
						for (l = 0; l < 8; l++) {
							if (bP[j] == octreeMesh.e[cE2[j][k]][l]) {
								testP[l][0] = target[0]; testP[l][1] = target[1]; testP[l][2] = target[2];
							}
							else {
								testP[l][0] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][0];
								testP[l][1] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][1];
								testP[l][2] = octreeMesh.v[octreeMesh.e[cE2[j][k]][l]][2];
							}
						}
						if (Sj(testP[0], testP[1], testP[2], testP[3], testP[4], testP[5], testP[6], testP[7]) <= ELEM_THRES) {
							pair = false;
							break;
						}
					}
					if (pair) {
						octreeMesh.v[bP[j]][0] = target[0];
						octreeMesh.v[bP[j]][1] = target[1];
						octreeMesh.v[bP[j]][2] = target[2];
					}
				}

				for (j = 0; j < sPIdx; j++) {
					tmp[0] = 0; tmp[1] = 0; tmp[2] = 0;
					for (k = 0; k < cP[j].size(); k++) {
						tmp[0] += octreeMesh.v[cP[j][k]][0];
						tmp[1] += octreeMesh.v[cP[j][k]][1];
						tmp[2] += octreeMesh.v[cP[j][k]][2];
					}
					tmp[0] /= cP[j].size(); tmp[1] /= cP[j].size(); tmp[2] /= cP[j].size();
					prop = rand() / (RAND_MAX * 1.0f);
					tmp[0] = prop * tmp[0] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][0];
					tmp[1] = prop * tmp[1] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][1];
					tmp[2] = prop * tmp[2] + (1 - prop) * octreeMesh.v[bP[j + sPIdx]][2];

					minDist = MAX_NUM2;
					for (k = 0; k < triMesh.eNum; k++) {
						dis = PointToTri(triMesh.v[triMesh.e[k][0]], triMesh.v[triMesh.e[k][1]], triMesh.v[triMesh.e[k][2]],
							tmp, tmp2, minDist);
						if (dis < minDist) {
							minDist = dis;
							target[0] = tmp2[0]; target[1] = tmp2[1]; target[2] = tmp2[2];
						}
					}

					pair = true;
					for (k = 0; k < cE[j].size(); k++) {
						for (l = 0; l < 8; l++) {
							if (bP[j + sPIdx] == octreeMesh.e[cE[j][k]][l]) {
								testP[l][0] = target[0]; testP[l][1] = target[1]; testP[l][2] = target[2];
							}
							else {
								testP[l][0] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][0];
								testP[l][1] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][1];
								testP[l][2] = octreeMesh.v[octreeMesh.e[cE[j][k]][l]][2];
							}
						}
						if (Sj(testP[0], testP[1], testP[2], testP[3], testP[4], testP[5], testP[6], testP[7]) <= ELEM_THRES) {
							pair = false;
							break;
						}
					}
					if (pair) {
						octreeMesh.v[bP[j + sPIdx]][0] = target[0];
						octreeMesh.v[bP[j + sPIdx]][1] = target[1];
						octreeMesh.v[bP[j + sPIdx]][2] = target[2];
					}
				}
				octreeMesh.WriteToVtk("finalMesh.vtk", BOX_LENGTH_RATIO, START_POINT);
				smallDist = 114514;
			}
		}
	}
}

inline void hexGen::InitiateElementValence() {
	// face index: bottom, front, left, right, back, top
	// points in each face are arranged: [small, small], [big, small], [big, big], [small, big]
	// 0123, 0154, 0374, 1265, 3267, 4567

	// initialization
	int i, j, k;

	elementValence = new int[leafNum][6][4];
	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 6; j++)
			for (k = 0; k < 4; k++)
				elementValence[i][j][k] = -1;

	elementValenceNumber = new int[leafNum][6];
	for (i = 0; i < leafNum; i++)
		for (j = 0; j < 6; j++)
			elementValenceNumber[i][j] = 0;

	int fIdTmp[2][4];// 0: checking face, 1: pairing face
	int l, m, pFlag;
	// pFlag: -1 means no pairing; 0-3 means the pairing index; 4 means all paired  
	// pair faces
	for (i = 0; i < leafNum - 1; i++)
		for (j = 0; j < 6; j++) {
			for (k = 0; k < 4; k++)// checking face index
				fIdTmp[0][k] = octreeMesh.e[i][fId[j][k]];
			for (k = i + 1; k < leafNum; k++) {
				pFlag = -1;
				l = 5 - j;// the opposite face
				for (m = 0; m < 4; m++) {// pairing face index
					fIdTmp[1][m] = octreeMesh.e[k][fId[l][m]];
					if (fIdTmp[0][m] == fIdTmp[1][m])// find a paired point
						if (pFlag == -1)
							pFlag = m;
						else {
							pFlag = 4;
							break;
						}
				}
				if (pFlag == 4) {// the checking face and the pairing face have the same size
					elementValence[i][j][0] = k;
					elementValence[k][l][0] = i;

					elementValenceNumber[i][j] = 1;
					elementValenceNumber[k][l] = 1;
				}
				else if (pFlag != -1) {// the checking face and the pairing face have different face size
					elementValence[i][j][pFlag] = k;
					elementValence[k][l][pFlag] = i;
					if (elementValenceNumber[i][j] == 0)
						elementValence[i][j][0] = k;
					if (elementValenceNumber[k][l] == 0)
						elementValence[k][l][0] = i;
					elementValenceNumber[i][j]++;
					elementValenceNumber[k][l]++;
				}
			}
		}
}
