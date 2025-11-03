/* ========================================================================
 * inpolygon_mex_v23.cpp
 *
 * This version FIXES the catastrophic algorithmic failure of all
 * previous C++ scan-line versions (v16-v22).
 *
 * (FIX 24) The previous versions incorrectly incremented diffCount[0]
 * AND flipped the 'isInside' parity bit for intersections
 * to the left of the grid (OOB-Left), causing a double-count
 * and totally incorrect results.
 *
 * This version introduces a critical 'else if' block to
 * strictly separate the logic:
 * 1. OOB-Left intersections ONLY affect 'isInside' (parity).
 * 2. In-Grid intersections ONLY affect 'diffCount' (cumsum).
 *
 * This is the correct, robust, and fast scan-line implementation.
 *
 * COMPILATION (Windows, with OpenMP):
 * mex -v CXXFLAGS="$CXXFLAGS /std:c++17 /O2 /openmp" LINKLIBS="$LINKLIBS /openmp" inpolygon_mex_v23.cpp
 *
 * COMPILATION (Linux/macOS, with OpenMP):
 * mex -v CXXFLAGS="$CXXFLAGS -std:c++17 -O2 -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" inpolygon_mex_v23.cpp
 *
 * ========================================================================*/

#include "mex.h"
#include <vector>
#include <cmath> // For std::isnan, std::isinf
#include <algorithm>
#include <limits>

// --- OpenMP include ---
#ifdef _OPENMP
#include <omp.h>
#endif

// --- Robust floating-point constants ---
constexpr double EPSILON = 1e-9;
constexpr double EPSILON_SQ = 1e-18;

// --- Helper Structure for Points ---
struct Point {
    double x, y;
};

// --- Helper: Robust On-Segment Check ---
bool isPointOnSegment(const Point& p, const Point& a, const Point& b) {
    double edge_len_sq = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
    if (edge_len_sq < EPSILON_SQ) {
        return (std::abs(p.x - a.x) < EPSILON && std::abs(p.y - a.y) < EPSILON);
    }
    
    double cross_product = (p.y - a.y) * (b.x - a.x) - (p.x - a.x) * (b.y - a.y);
    double distance_sq = (cross_product * cross_product) / edge_len_sq;

    if (distance_sq > EPSILON_SQ) {
        return false; 
    }

    double dot_product = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
    return (dot_product >= -EPSILON_SQ) && (dot_product <= edge_len_sq + EPSILON_SQ);
}

// --- Helper: Find grid indices for a world-coordinate BBox ---
void findSubGrid(const double* v_coord, size_t v_len, double v_min, double v_max, 
                 size_t& v_start_idx, size_t& v_end_idx) 
{
    const double* start_ptr = std::lower_bound(v_coord, v_coord + v_len, v_min - EPSILON);
    v_start_idx = (start_ptr - v_coord);

    const double* end_ptr = std::upper_bound(v_coord + v_start_idx, v_coord + v_len, v_max + EPSILON);
    v_end_idx = (end_ptr - v_coord); 
}

// --- Helper: Find the *next* grid index (for diffCount) ---
size_t findNextGridIndex(const double* v_coord, size_t v_len, double intersect_val) {
    // std::upper_bound finds the first element *greater* than intersect_val
    const double* ptr = std::upper_bound(v_coord, v_coord + v_len, intersect_val - EPSILON);
    return (ptr - v_coord); // Returns index [0, v_len]
}


// --- Main MEX Function ---
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // --- 1. Check Inputs ---
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_mex:nrhs", "Four inputs required: value_x, value_y, xv, yv.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_mex:nlhs", "Two outputs required: [inMesh, onedgeMesh].");
    }
    if (mxGetM(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_mex:badInput1", "Input 1 (value_x) must be a (1, N) real row vector.");
    }
    if (mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_mex:badInput2", "Input 2 (value_y) must be an (M, 1) real column vector.");
    }

    // --- 2. Get Grid Vectors ---
    const double* value_x = mxGetPr(prhs[0]);
    const double* value_y = mxGetPr(prhs[1]);
    const size_t N = mxGetN(prhs[0]); // Columns (width)
    const size_t M = mxGetM(prhs[1]); // Rows (height)
    const size_t numGridPoints = M * N;

    // --- 3. Get Polygon Vectors ---
    const double* xv_full = mxGetPr(prhs[2]);
    const double* yv_full = mxGetPr(prhs[3]);
    const size_t V_full = mxGetNumberOfElements(prhs[2]);

    // --- 4. Allocate Outputs ---
    plhs[0] = mxCreateLogicalMatrix(M, N); // inMesh
    plhs[1] = mxCreateLogicalMatrix(M, N); // onedgeMesh
    bool* inMeshPtr = mxGetLogicals(plhs[0]);
    bool* onedgeMeshPtr = mxGetLogicals(plhs[1]);
    
    std::fill(inMeshPtr, inMeshPtr + numGridPoints, false);
    std::fill(onedgeMeshPtr, onedgeMeshPtr + numGridPoints, false);

    // --- 5. Allocate Intermediates ---
    std::vector<int> diffCount(numGridPoints, 0); 
    std::vector<bool> isInsideY(M, false); 
    std::vector<bool> isInsideX(N, false); 

    // --- 6. Find Polygon BBox and Scan Direction ---
    double poly_min_x = std::numeric_limits<double>::infinity();
    double poly_max_x = -std::numeric_limits<double>::infinity();
    double poly_min_y = std::numeric_limits<double>::infinity();
    double poly_max_y = -std::numeric_limits<double>::infinity();
    
    for(size_t i = 0; i < V_full; ++i) {
        if (!mxIsNaN(xv_full[i])) {
            poly_min_x = std::min(poly_min_x, xv_full[i]);
            poly_max_x = std::max(poly_max_x, xv_full[i]);
        }
        if (!mxIsNaN(yv_full[i])) {
            poly_min_y = std::min(poly_min_y, yv_full[i]);
            poly_max_y = std::max(poly_max_y, yv_full[i]);
        }
    }
    
    size_t r_start_bbox, r_end_bbox, c_start_bbox, c_end_bbox;
    findSubGrid(value_y, M, poly_min_y, poly_max_y, r_start_bbox, r_end_bbox);
    findSubGrid(value_x, N, poly_min_x, poly_max_x, c_start_bbox, c_end_bbox);
    
    size_t M_sub = (r_end_bbox > r_start_bbox) ? (r_end_bbox - r_start_bbox) : 0;
    size_t N_sub = (c_end_bbox > c_start_bbox) ? (c_end_bbox - c_start_bbox) : 0;

    bool useYScan = (M_sub <= N_sub); 
    
    // --- 7. Process Polygons (Handle NaNs) ---
    std::vector<size_t> nan_indices;
    nan_indices.push_back(static_cast<size_t>(-1));
    for (size_t i = 0; i < V_full; ++i) {
        if (mxIsNaN(xv_full[i]) || mxIsNaN(yv_full[i])) {
            nan_indices.push_back(i);
        }
    }
    nan_indices.push_back(V_full); 

    std::vector<Point> poly_pts; 
    
    for (size_t p = 0; p < nan_indices.size() - 1; ++p) {
        size_t start_idx = nan_indices[p] + 1;
        size_t end_idx = nan_indices[p+1] - 1;
        
        if (end_idx < start_idx) continue;
        
        poly_pts.clear();
        for(size_t i = start_idx; i <= end_idx; ++i) {
            poly_pts.push_back({xv_full[i], yv_full[i]});
        }
        
        if (poly_pts.front().x != poly_pts.back().x || poly_pts.front().y != poly_pts.back().y) {
            poly_pts.push_back(poly_pts.front());
        }
        
        size_t V = poly_pts.size();
        if (V < 3) continue;
        
        // --- 7.A: Part A (On-Edge) C++ Implementation ---
        #pragma omp parallel for
        for (size_t e = 0; e < V - 1; ++e) {
            Point a = poly_pts[e];
            Point b = poly_pts[e+1];
            
            double edge_len_sq = (b.x - a.x)*(b.x - a.x) + (b.y - a.y)*(b.y - a.y);
            if (edge_len_sq < EPSILON_SQ) continue;
            
            double minX = std::min(a.x, b.x) - EPSILON;
            double maxX = std::max(a.x, b.x) + EPSILON;
            double minY = std::min(a.y, b.y) - EPSILON;
            double maxY = std::max(a.y, b.y) + EPSILON;
            
            size_t r_start, r_end, c_start, c_end;
            findSubGrid(value_y, M, minY, maxY, r_start, r_end);
            findSubGrid(value_x, N, minX, maxX, c_start, c_end);

            for (size_t r = r_start; r < r_end; ++r) {
                for (size_t c = c_start; c < c_end; ++c) {
                    size_t idx = c * M + r; 
                    if (onedgeMeshPtr[idx]) continue; 
                    
                    if (isPointOnSegment({value_x[c], value_y[r]}, a, b)) {
                        #pragma omp atomic
                        onedgeMeshPtr[idx] = true;
                    }
                }
            }
        } // End Part A edge loop

        // --- 7.B: Part B (Scan-Line) C++ Implementation ---
        if (useYScan) {
            // --- Y-SCAN ---
            for (size_t e = 0; e < V - 1; ++e) {
                Point p1 = poly_pts[e];
                Point p2 = poly_pts[e+1];
                
                if (std::abs(p1.y - p2.y) < EPSILON) continue; // Skip horizontal
                
                double inv_slope = (p2.x - p1.x) / (p2.y - p1.y);
                
                double min_y = std::min(p1.y, p2.y);
                double max_y = std::max(p1.y, p2.y);
                size_t r_start, r_end;
                findSubGrid(value_y, M, min_y, max_y, r_start, r_end);

                for (size_t r = r_start; r < r_end; ++r) {
                    double y = value_y[r];
                    
                    if (y < min_y - EPSILON || y >= max_y - EPSILON) continue;
                    
                    double x_intersect = p1.x + (y - p1.y) * inv_slope;
                    
                    if (std::isnan(x_intersect) || std::isinf(x_intersect)) {
                        continue;
                    }
                    
                    // --- FIX 24: Use 'else if' to separate OOB-Left from In-Grid ---
                    if (x_intersect < value_x[0] - EPSILON) {
                        // OOB-Left: Only flip parity
                        isInsideY[r] = !isInsideY[r];
                    } 
                    else if (x_intersect <= value_x[N-1] + EPSILON) {
                        // In-Grid: Find *next* index and increment diffCount
                        size_t c = findNextGridIndex(value_x, N, x_intersect);
                        if (c < N) {
                            diffCount[c * M + r]++; 
                        }
                    }
                    // OOB-Right (x_intersect > value_x[N-1]) does nothing.
                }
            }
        } else {
            // --- X-SCAN ---
            for (size_t e = 0; e < V - 1; ++e) {
                Point p1 = poly_pts[e];
                Point p2 = poly_pts[e+1];
                
                if (std::abs(p1.x - p2.x) < EPSILON) continue; // Skip vertical
                
                double slope = (p2.y - p1.y) / (p2.x - p1.x);

                double min_x = std::min(p1.x, p2.x);
                double max_x = std::max(p1.x, p2.x);
                size_t c_start, c_end;
                findSubGrid(value_x, N, min_x, max_x, c_start, c_end);
                
                for (size_t c = c_start; c < c_end; ++c) {
                    double x = value_x[c];
                    
                    if (x < min_x - EPSILON || x >= max_x - EPSILON) continue;
                    
                    double y_intersect = p1.y + (x - p1.x) * slope;
                    
                    if (std::isnan(y_intersect) || std::isinf(y_intersect)) {
                        continue;
                    }
                    
                    // --- FIX 24: Use 'else if' to separate OOB-Below from In-Grid ---
                    if (y_intersect < value_y[0] - EPSILON) {
                        // OOB-Below: Only flip parity
                        isInsideX[c] = !isInsideX[c];
                    } 
                    else if (y_intersect <= value_y[M-1] + EPSILON) {
                        // In-Grid: Find *next* index and increment diffCount
                        size_t r = findNextGridIndex(value_y, M, y_intersect);
                        if (r < M) {
                            diffCount[c * M + r]++;
                        }
                    }
                    // OOB-Above (y_intersect > value_y[M-1]) does nothing.
                }
            }
        } // End Part B scan-direction
    } // End polygon loop
    
    // --- 8. Finalize Outputs (Now parallelized!) ---
    if (useYScan) {
        // Y-Scan: cumsum along dimension 2 (X)
        #pragma omp parallel for
        for (size_t r = 0; r < M; ++r) {
            int crossingCount = 0;
            bool parity = isInsideY[r]; 
            for (size_t c = 0; c < N; ++c) {
                size_t idx = c * M + r; 
                crossingCount += diffCount[idx];
                bool totalInclusion = (crossingCount % 2) != 0;
                totalInclusion = totalInclusion ^ parity;
                
                if (!onedgeMeshPtr[idx]) { 
                    inMeshPtr[idx] = totalInclusion;
                }
            }
        }
    } else {
        // X-Scan: cumsum along dimension 1 (Y)
        #pragma omp parallel for
        for (size_t c = 0; c < N; ++c) {
            int crossingCount = 0;
            bool parity = isInsideX[c]; 
            for (size_t r = 0; r < M; ++r) {
                size_t idx = c * M + r;
                crossingCount += diffCount[idx];
                bool totalInclusion = (crossingCount % 2) != 0;
                totalInclusion = totalInclusion ^ parity;
                
                if (!onedgeMeshPtr[idx]) {
                    inMeshPtr[idx] = totalInclusion;
                }
            }
        }
    }
}