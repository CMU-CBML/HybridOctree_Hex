#ifndef INITIALIZATION_H
#define INITIALIZATION_H

// constants
const double PI = 3.1415926535897932384626433;
const int VOXEL_SIZE = 10;// log2 voxel size
const double C_THRES[5] = { 0, 0, 0.4, 0.8, 1.6};// 4 5 6 7 8
const double H_THRES[5] = { 16, 8, 4, 2, 1};// 4 5 6 7 8
const int MAX_NUM = 100000000;// maximum number
const int MAX_NUM2 = 2147483647;
const double DIST_THRES = 1e-12;// threshold when judging point overlap
const double DIST_THRES2 = sqrt(DIST_THRES);
const double DIST_THRES3 = 1.0f / (1 << (2 * VOXEL_SIZE - 6));
const double OUT_IN_RATIO = 0.15;
const double COS_THRES = 1e-2;// eliminate invalid unmanifold hexes, set it slightly greater than zero to avoid numerical error problem
const int UPDATE_EVERY = 1000;
const int SMOOTH_EPOCH = 1;// smoothing optimization epoch
const double CELL_DETECT = 1.0;

// defines
#define ORIENT_2D(a, b, c)  ((a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]))

#define CROSS(dest, v1, v2)   dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
                              dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
                              dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2)   dest[0] = v1[0] - v2[0]; \
                            dest[1] = v1[1] - v2[1]; \
                            dest[2] = v1[2] - v2[2];

#define SCALAR(dest, alpha, v)   dest[0] = alpha * v[0]; \
                                 dest[1] = alpha * v[1]; \
                                 dest[2] = alpha * v[2];

#define CHECK_MIN_MAX(p1, q1, r1, p2, q2, r2) \
{ \
    SUB(v1, p2, q1) \
    SUB(v2, p1, q1) \
    CROSS(N1, v1, v2) \
    SUB(v1, q2, q1) \
    if (DOT(v1, N1) >= -DIST_THRES) return 0; \
    SUB(v1, p2, p1) \
    SUB(v2, r1, p1) \
    CROSS(N1, v1, v2) \
    SUB(v1, r2, p1) \
    if (DOT(v1, N1) >= -DIST_THRES) return 0; \
    return 1; \
}

// permutation in a canonical form of T2's vertices
#define TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2) \
{ \
    if (dp2 > 0) { \
        if      (dq2 > 0) CHECK_MIN_MAX(p1, r1, q1, r2, p2, q2) \
        else if (dr2 > 0) CHECK_MIN_MAX(p1, r1, q1, q2, r2, p2) \
        else              CHECK_MIN_MAX(p1, q1, r1, p2, q2, r2) \
    } else if (dp2 < 0) { \
        if      (dq2 < 0) CHECK_MIN_MAX(p1, q1, r1, r2, p2, q2) \
        else if (dr2 < 0) CHECK_MIN_MAX(p1, q1, r1, q2, r2, p2) \
        else              CHECK_MIN_MAX(p1, r1, q1, p2, q2, r2) \
    } else { \
        if (dq2 < 0) { \
            if (dr2 >= 0) CHECK_MIN_MAX(p1, r1, q1, q2, r2, p2) \
            else          CHECK_MIN_MAX(p1, q1, r1, p2, q2, r2) \
        } else if (dq2 > 0) { \
            if (dr2 > 0) CHECK_MIN_MAX(p1, r1, q1, p2, q2, r2) \
            else         CHECK_MIN_MAX(p1, q1, r1, q2, r2, p2) \
        } else { \
            if      (dr2 > 0) CHECK_MIN_MAX(p1, q1, r1, r2, p2, q2) \
            else if (dr2 < 0) CHECK_MIN_MAX(p1, r1, q1, r2, p2, q2) \
            else return Coplanar(p1, q1, r1, p2, q2, r2, N1, N2); \
        } \
    } \
}

#define CONSTRUCT_INTERSECTION(p1, q1, r1, p2, q2, r2) \
{ \
    SUB(v1, q1, p1) \
    SUB(v2, r2, p1) \
    CROSS(N, v1, v2) \
    SUB(v, p2, p1) \
    if (DOT(v, N) > 0) { \
        SUB(v1, r1, p1) \
        CROSS(N, v1, v2) \
        if (DOT(v, N) <= 0) { \
            SUB(v2, q2, p1) \
            CROSS(N, v1, v2) \
            if (DOT(v, N) > 0) { \
                SUB(v1, p1, p2) \
                SUB(v2, p1, r1) \
                alpha = DOT(v1, N2) / DOT(v2, N2); \
                SCALAR(v1, alpha, v2) \
                SUB(source, p1, v1) \
                SUB(v1, p2, p1) \
                SUB(v2, p2, r2) \
                alpha = DOT(v1, N1) / DOT(v2, N1); \
                SCALAR(v1, alpha, v2) \
                SUB(target, p2, v1) \
                return 1; \
            } else { \
                SUB(v1, p2, p1) \
                SUB(v2, p2, q2) \
                alpha = DOT(v1, N1) / DOT(v2, N1); \
                SCALAR(v1, alpha, v2) \
                SUB(source, p2, v1) \
                SUB(v1, p2, p1) \
                SUB(v2, p2, r2) \
                alpha = DOT(v1, N1) / DOT(v2, N1); \
                SCALAR(v1, alpha, v2) \
                SUB(target, p2, v1) \
                return 1; \
            } \
        } else { \
            return 0; \
        } \
    } else { \
        SUB(v2, q2, p1) \
        CROSS(N, v1, v2) \
        if (DOT(v, N) < 0) { \
            return 0; \
        } else { \
            SUB(v1, r1, p1) \
            CROSS(N, v1, v2) \
            if (DOT(v, N) >= 0) { \
                SUB(v1, p1, p2) \
                SUB(v2, p1, r1) \
                alpha = DOT(v1, N2) / DOT(v2, N2); \
                SCALAR(v1, alpha, v2) \
                SUB(source, p1, v1) \
                SUB(v1, p1, p2) \
                SUB(v2, p1, q1) \
                alpha = DOT(v1, N2) / DOT(v2, N2); \
                SCALAR(v1, alpha, v2) \
                SUB(target, p1, v1) \
                return 1; \
            } else { \
                SUB(v1, p2, p1) \
                SUB(v2, p2, q2) \
                alpha = DOT(v1, N1) / DOT(v2, N1); \
                SCALAR(v1, alpha, v2) \
                SUB(source, p2, v1) \
                SUB(v1, p1, p2) \
                SUB(v2, p1, q1) \
                alpha = DOT(v1, N2) / DOT(v2, N2); \
                SCALAR(v1, alpha, v2) \
                SUB(target, p1, v1) \
                return 1; \
            } \
        } \
    } \
}

#define TRI_TRI_INTER_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2) \
{ \
    if (dp2 > 0) { \
        if      (dq2 > 0) CONSTRUCT_INTERSECTION(p1, r1, q1, r2, p2, q2) \
        else if (dr2 > 0) CONSTRUCT_INTERSECTION(p1, r1, q1, q2, r2, p2) \
        else              CONSTRUCT_INTERSECTION(p1, q1, r1, p2, q2, r2) \
    } else if (dp2 < 0) { \
        if      (dq2 < 0) CONSTRUCT_INTERSECTION(p1, q1, r1, r2, p2, q2) \
        else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1, q1, r1, q2, r2, p2) \
        else              CONSTRUCT_INTERSECTION(p1, r1, q1, p2, q2, r2) \
    } else { \
        if (dq2 < 0) { \
            if (dr2 >= 0) CONSTRUCT_INTERSECTION(p1, r1, q1, q2, r2, p2) \
            else          CONSTRUCT_INTERSECTION(p1, q1, r1, p2, q2, r2) \
        } else if (dq2 > 0) { \
            if (dr2 > 0) CONSTRUCT_INTERSECTION(p1, r1, q1, p2, q2, r2) \
            else         CONSTRUCT_INTERSECTION(p1, q1, r1, q2, r2, p2) \
        } else  { \
            if      (dr2 > 0) CONSTRUCT_INTERSECTION(p1, q1, r1, r2, p2, q2) \
            else if (dr2 < 0) CONSTRUCT_INTERSECTION(p1, r1, q1, r2, p2, q2) \
            else { \
                *coplanar = 1; \
                return Coplanar(p1, q1, r1, p2, q2, r2, N1, N2);\
            } \
        } \
    } \
}

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) \
{ \
    if (ORIENT_2D(R2, P2, Q1) >= 0) { \
        if (ORIENT_2D(R2, Q2, Q1) <= 0) { \
            if (ORIENT_2D(P1, P2, Q1) > 0) { \
                if (ORIENT_2D(P1, Q2, Q1) <= 0) return 1; \
            } else { \
                if (ORIENT_2D(P1, P2, R1) >= 0) { \
                    if (ORIENT_2D(Q1, R1, P2) >= 0) return 1; \
                } \
            } \
        } else { \
            if (ORIENT_2D(P1, Q2, Q1) <= 0) { \
                if (ORIENT_2D(R2, Q2, R1) <= 0) { \
                    if (ORIENT_2D(Q1, R1, Q2) >= 0) return 1; \
                } \
            } \
        } \
    } else { \
        if (ORIENT_2D(R2, P2, R1) >= 0) { \
            if (ORIENT_2D(Q1, R1, R2) >= 0) { \
                if (ORIENT_2D(P1, P2, R1) >= 0) return 1;\
            } else { \
                if (ORIENT_2D(Q1, R1, Q2) >= 0) { \
                    if (ORIENT_2D(R2, R1, Q2) >= 0) return 1; \
                } \
            } \
        } \
    } \
    return 0; \
}

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) \
{ \
    if (ORIENT_2D(R2, P2, Q1) >= 0) { \
        if (ORIENT_2D(P1, P2, Q1) >= 0) { \
            if (ORIENT_2D(P1, Q1, R2) >= 0) return 1; \
        } else { \
            if (ORIENT_2D(Q1, R1, P2) >= 0) { \
                if (ORIENT_2D(R1, P1, P2) >= 0) return 1; \
            } \
        } \
    } else { \
        if (ORIENT_2D(R2, P2, R1) >= 0) { \
            if (ORIENT_2D(P1, P2, R1) >= 0) { \
                if (ORIENT_2D(P1, R1, R2) >= 0) return 1; \
                if (ORIENT_2D(Q1, R1, R2) >= 0) return 1; \
            } \
        } \
    } \
    return 0; \
}

#endif INITIALIZATION_H