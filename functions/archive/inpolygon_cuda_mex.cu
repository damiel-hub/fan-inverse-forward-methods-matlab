/* ========================================================================
 * inpolygon_cuda_mex.cu
 *
 * This is a CUDA MEX file for a massively parallel point-in-polygon
 * check. It uses the ray-casting algorithm, which is GPU-friendly.
 *
 * This implementation is a rewrite of the CPU-based scan-line algorithm
 * and is designed to run on an NVIDIA GPU.
 *
 * COMPILATION (Windows/Linux/macOS):
 * 1. Make sure you have the NVIDIA CUDA Toolkit installed.
 * 2. In MATLAB, run:
 * mexcuda inpolygon_cuda_mex.cu
 *
 * ========================================================================*/

#include "mex.h"
#include "cuda_runtime.h"
#include <cmath>
#include <limits>
#include <stdexcept> // For exception handling
#include <vector>    // <--- ADD THIS INCLUDE

// --- CUDA Error-Checking Macro ---
// This macro will check for CUDA errors and report them to the MATLAB console.
static void handleCudaError(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:cudaError",
                          "CUDA error in %s at line %d: %s",
                          file, line, cudaGetErrorString(err));
    }
}
#define CUDA_CHECK(err) (handleCudaError(err, __FILE__, __LINE__))

// --- Robust floating-point constants (for GPU) ---
__device__ constexpr double EPSILON = 1e-9;
__device__ constexpr double EPSILON_SQ = 1e-18;

// --- Helper: Robust On-Segment Check (GPU Device Function) ---
// This function runs on the GPU.
__device__ bool isPointOnSegment_device(
    double px, double py, 
    double ax, double ay, 
    double bx, double by) 
{
    double edge_len_sq = (bx - ax)*(bx - ax) + (by - ay)*(by - ay);
    
    // Check for a zero-length edge (just a point)
    if (edge_len_sq < EPSILON_SQ) {
        return (std::abs(px - ax) < EPSILON && std::abs(py - ay) < EPSILON);
    }
    
    // --- 1. Collinearity Check ---
    // Calculate the cross product. If it's non-zero, the point is not on the line.
    double cross_product = (py - ay) * (bx - ax) - (px - ax) * (by - ay);
    
    // Use squared distance for efficiency (avoids sqrt)
    double distance_sq = (cross_product * cross_product) / edge_len_sq;

    if (distance_sq > EPSILON_SQ) {
        return false; // Point is not collinear
    }

    // --- 2. Bounding Box Check ---
    // Check if the point is within the segment's "reach" using dot product.
    double dot_product = (px - ax) * (bx - ax) + (py - ay) * (by - ay);
    
    // Point is on the line, but is it on the segment?
    // It is if 0 <= dot_product <= edge_len_sq
    return (dot_product >= -EPSILON_SQ) && (dot_product <= edge_len_sq + EPSILON_SQ);
}


// --- CUDA Kernel: Ray Casting (MODIFIED FOR MULTI-POLYGON) ---
// This function is the core logic that runs on the GPU.
// Each thread executes this function for one grid point.
__global__ void rayCastKernel(
    const double* value_x,  // [N]
    const double* value_y,  // [M]
    const double* xv,       // [V_clean_total] (Concatenated vertices)
    const double* yv,       // [V_clean_total]
    bool* inMesh,           // [M x N]
    bool* onMesh,           // [M x N]
    const size_t* poly_offsets, // [numPolygons + 1]
    size_t M,               // Grid height
    size_t N,               // Grid width
    size_t numPolygons)     // Number of polygons
{
    // --- 1. Get this thread's unique ID ---
    size_t idx = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    size_t numGridPoints = M * N;

    // --- 2. Check if thread is out of bounds ---
    if (idx >= numGridPoints) {
        return;
    }

    // --- 3. Convert 1D index (idx) to 2D grid coordinates (r, c) ---
    size_t r = idx % M;
    size_t c = idx / M;
    double px = value_x[c];
    double py = value_y[r];

    // --- 4. Initialize thread-local variables ---
    int totalCrossingCount = 0;
    bool totalOnEdge = false;

    // --- 5. Loop over all polygons ---
    for (size_t p = 0; p < numPolygons; ++p) {
        
        // Find the start and end vertex indices for *this* polygon
        size_t poly_start_idx = poly_offsets[p];
        size_t poly_end_idx   = poly_offsets[p+1];
        size_t numVerticesInPoly = poly_end_idx - poly_start_idx;

        // A polygon must have at least 3 vertices (e.g., v0, v1, v2) 
        // which means at least 3 edges (v0-v1, v1-v2, v2-v0).
        // Since our data is pre-closed, 3 vertices = 2 edges (v0-v1, v1-v0) -
        // let's require 3 edges, so 4 vertices.
        if (numVerticesInPoly < 4) { // e.g., v0, v1, v2, v0
            continue;
        }
        
        size_t numEdgesInPoly = numVerticesInPoly - 1; // It's pre-closed

        int polyCrossingCount = 0;
        bool polyOnEdge = false;

        // --- 5.A. Loop over all edges *in this polygon* ---
        for (size_t e_local = 0; e_local < numEdgesInPoly; ++e_local) {
            size_t e_global = poly_start_idx + e_local;
            
            double ax = xv[e_global];
            double ay = yv[e_global];
            double bx = xv[e_global + 1];
            double by = yv[e_global + 1];

            // --- 5.B. Check if the point is *on* this edge ---
            if (isPointOnSegment_device(px, py, ax, ay, bx, by)) {
                polyOnEdge = true;
                break; // No need to check other edges in *this* polygon
            }

            // --- 5.C. Ray-Casting Check ---
            bool y_span = ((ay <= py && py < by) || (by <= py && py < ay));
            
            if (y_span) {
                double x_intersect = (py - ay) * (bx - ax) / (by - ay) + ax;
                if (px < x_intersect) {
                    polyCrossingCount++;
                }
            }
        } // End edge loop

        // --- 5.D. Combine this polygon's results ---
        if (polyOnEdge) {
            totalOnEdge = true;
            totalCrossingCount = 0; // 'on' overrides 'in'
            break; // No need to check other polygons
        } else {
            // XOR logic is handled by sum and modulo 2 at the end
            totalCrossingCount += polyCrossingCount;
        }
    } // End polygon loop

    // --- 6. Write final results ---
    onMesh[idx] = totalOnEdge;
    
    if (totalOnEdge) {
        inMesh[idx] = false; // 'on' takes precedence
    } else {
        // Point is 'in' if it has an ODD total number of crossings.
        inMesh[idx] = (totalCrossingCount % 2) != 0;
    }
}


// --- Main MEX Function (CPU Host) ---
// This is the entry point that MATLAB calls.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // --- 1. Check Inputs ---
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:nrhs", "Four inputs required: value_x, value_y, xv, yv.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:nlhs", "Two outputs required: [inMesh, onedgeMesh].");
    }
    // (Add more input validation as in your original...)

    // --- 2. Get Grid Vectors (CPU) ---
    const double* value_x_cpu = mxGetPr(prhs[0]);
    const double* value_y_cpu = mxGetPr(prhs[1]);
    const size_t N = mxGetN(prhs[0]); // Columns (width)
    const size_t M = mxGetM(prhs[1]); // Rows (height)
    const size_t numGridPoints = M * N;

    // --- 3. Get Polygon Vectors (CPU) ---
    const double* xv_full_cpu = mxGetPr(prhs[2]);
    const double* yv_full_cpu = mxGetPr(prhs[3]);
    const size_t V_full = mxGetNumberOfElements(prhs[2]);

    if (V_full != mxGetNumberOfElements(prhs[3])) {
         mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:xyMismatch", "xv and yv must have the same number of elements.");
    }
    
    // --- 3.B. Pre-process Polygons (Handle NaNs) ---
    // This logic is adapted from your original CPU v23 code.
    std::vector<size_t> nan_indices;
    std::vector<double> xv_clean_vec;
    std::vector<double> yv_clean_vec;
    std::vector<size_t> poly_offsets_vec;

    nan_indices.push_back(static_cast<size_t>(-1)); // Start before the first element
    for (size_t i = 0; i < V_full; ++i) {
        if (mxIsNaN(xv_full_cpu[i]) || mxIsNaN(yv_full_cpu[i])) {
            nan_indices.push_back(i);
        }
    }
    nan_indices.push_back(V_full); // End after the last element
    
    for (size_t p = 0; p < nan_indices.size() - 1; ++p) {
        size_t start_idx = nan_indices[p] + 1;
        size_t end_idx = nan_indices[p+1] - 1;
        
        // Check for empty polygon parts (e.g., consecutive NaNs)
        if (end_idx < start_idx || (end_idx - start_idx + 1) < 2) continue;
        
        // Store the starting index *in the clean vector*
        poly_offsets_vec.push_back(xv_clean_vec.size());
        
        // Add all vertices for this polygon part
        for(size_t i = start_idx; i <= end_idx; ++i) {
            xv_clean_vec.push_back(xv_full_cpu[i]);
            yv_clean_vec.push_back(yv_full_cpu[i]);
        }
        
        // Check if the polygon is closed. If not, close it.
        double first_x = xv_clean_vec[poly_offsets_vec.back()];
        double first_y = yv_clean_vec[poly_offsets_vec.back()];
        double last_x  = xv_clean_vec.back();
        double last_y  = yv_clean_vec.back();

        if (std::abs(first_x - last_x) > EPSILON || std::abs(first_y - last_y) > EPSILON) {
            xv_clean_vec.push_back(first_x);
            yv_clean_vec.push_back(first_y);
        }
    }
    // Add the final offset, which is the total number of clean vertices
    poly_offsets_vec.push_back(xv_clean_vec.size());
    
    const size_t numPolygons = poly_offsets_vec.size() - 1;
    const size_t V_clean_total = xv_clean_vec.size();

    if (numPolygons == 0 || V_clean_total == 0) {
        // No valid polygons were found, just return empty logicals
        plhs[0] = mxCreateLogicalMatrix(M, N);
        plhs[1] = mxCreateLogicalMatrix(M, N);
        return;
    }


    // --- 4. Allocate Outputs (CPU) ---
    plhs[0] = mxCreateLogicalMatrix(M, N); // inMesh
    plhs[1] = mxCreateLogicalMatrix(M, N); // onedgeMesh
    bool* inMesh_cpu = mxGetLogicals(plhs[0]);
    bool* onedgeMesh_cpu = mxGetLogicals(plhs[1]);

    // --- 5. Allocate Memory on GPU ---
    double* value_x_gpu;
    double* value_y_gpu;
    double* xv_gpu;
    double* yv_gpu;
    size_t* poly_offsets_gpu; // <-- ADD THIS
    bool* inMesh_gpu;
    bool* onMesh_gpu;

    CUDA_CHECK(cudaMalloc(&value_x_gpu, N * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&value_y_gpu, M * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&xv_gpu, V_clean_total * sizeof(double))); // <-- Use V_clean_total
    CUDA_CHECK(cudaMalloc(&yv_gpu, V_clean_total * sizeof(double))); // <-- Use V_clean_total
    CUDA_CHECK(cudaMalloc(&poly_offsets_gpu, (numPolygons + 1) * sizeof(size_t))); // <-- ALLOCATE
    CUDA_CHECK(cudaMalloc(&inMesh_gpu, numGridPoints * sizeof(bool)));
    CUDA_CHECK(cudaMalloc(&onMesh_gpu, numGridPoints * sizeof(bool)));

    // --- 6. Copy Data from CPU (Host) to GPU (Device) ---
    CUDA_CHECK(cudaMemcpy(value_x_gpu, value_x_cpu, N * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(value_y_gpu, value_y_cpu, M * sizeof(double), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(xv_gpu, xv_clean_vec.data(), V_clean_total * sizeof(double), cudaMemcpyHostToDevice)); // <-- Use clean vec
    CUDA_CHECK(cudaMemcpy(yv_gpu, yv_clean_vec.data(), V_clean_total * sizeof(double), cudaMemcpyHostToDevice)); // <-- Use clean vec
    CUDA_CHECK(cudaMemcpy(poly_offsets_gpu, poly_offsets_vec.data(), (numPolygons + 1) * sizeof(size_t), cudaMemcpyHostToDevice)); // <-- COPY
    // Initialize output arrays to false (optional, but good practice)
    CUDA_CHECK(cudaMemset(inMesh_gpu, 0, numGridPoints * sizeof(bool)));
    CUDA_CHECK(cudaMemset(onMesh_gpu, 0, numGridPoints * sizeof(bool)));


    // --- 7. Set up Kernel Launch Parameters ---
    // Threads per block (a common, reasonable size)
    int threadsPerBlock = 256; 
    
    // Number of blocks needed to cover all grid points
    // (We use integer division ceiling: (N + D - 1) / D)
    int numBlocks = (numGridPoints + threadsPerBlock - 1) / threadsPerBlock;

    // --- 8. Launch the CUDA Kernel ---
    rayCastKernel<<<numBlocks, threadsPerBlock>>>(
        value_x_gpu, value_y_gpu,
        xv_gpu, yv_gpu,
        inMesh_gpu, onMesh_gpu,
        poly_offsets_gpu,     // <-- PASS
        M, N,
        numPolygons           // <-- PASS
    );
    
    // Check for any errors during kernel launch
    CUDA_CHECK(cudaGetLastError());
    
    // Wait for the GPU to finish all work
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- 9. Copy Results from GPU (Device) to CPU (Host) ---
    CUDA_CHECK(cudaMemcpy(inMesh_cpu, inMesh_gpu, numGridPoints * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(onedgeMesh_cpu, onMesh_gpu, numGridPoints * sizeof(bool), cudaMemcpyDeviceToHost)); // <-- FIX: onedgeMesh_cpu

    // --- 10. Clean up GPU Memory ---
    CUDA_CHECK(cudaFree(value_x_gpu));
    CUDA_CHECK(cudaFree(value_y_gpu));
    CUDA_CHECK(cudaFree(xv_gpu));
    CUDA_CHECK(cudaFree(yv_gpu));
    CUDA_CHECK(cudaFree(poly_offsets_gpu)); // <-- FREE
    CUDA_CHECK(cudaFree(inMesh_gpu));
    CUDA_CHECK(cudaFree(onMesh_gpu));
}

