/* ========================================================================
 * inpolygon_cuda_mex_single.cu
 *
 * This is a CUDA MEX file for a massively parallel point-in-polygon
 * check. It uses the ray-casting algorithm, which is GPU-friendly.
 *
 * This version is modified to ONLY accept SINGLE PRECISION (float) inputs.
 *
 * COMPILATION (Windows/Linux/macOS):
 * 1. Make sure you have the NVIDIA CUDA Toolkit installed.
 * 2. In MATLAB, run:
 * mexcuda inpolygon_cuda_mex_single.cu
 *
 * ========================================================================*/
#include "mex.h"
#include "cuda_runtime.h"
#include <cmath>     // For std::isnan, std::abs
#include <limits>
#include <stdexcept> 
#include <vector>    

// --- CUDA Error-Checking Macro ---
static void handleCudaError(cudaError_t err, const char *file, int line) {
    if (err != cudaSuccess) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:cudaError",
                          "CUDA error in %s at line %d: %s",
                          file, line, cudaGetErrorString(err));
    }
}
#define CUDA_CHECK(err) (handleCudaError(err, __FILE__, __LINE__))

// --- Robust floating-point constants (for GPU) ---
__device__ constexpr float EPSILON_F = 1e-6f;
__device__ constexpr float EPSILON_SQ_F = 1e-12f; // (1e-6)^2

// --- Helper: Robust On-Segment Check (GPU Device Function) ---
// This function runs on the GPU.
__device__ bool isPointOnSegment_device(
    float px, float py, 
    float ax, float ay, 
    float bx, float by) 
{
    float edge_len_sq = (bx - ax)*(bx - ax) + (by - ay)*(by - ay);
    
    // Check for a zero-length edge (just a point)
    if (edge_len_sq < EPSILON_SQ_F) {
        return (fabsf(px - ax) < EPSILON_F && fabsf(py - ay) < EPSILON_F);
    }
    
    // --- 1. Collinearity Check ---
    float cross_product = (py - ay) * (bx - ax) - (px - ax) * (by - ay);
    
    float distance_sq = (cross_product * cross_product) / edge_len_sq;
    if (distance_sq > EPSILON_SQ_F) {
        return false; // Point is not collinear
    }
    // --- 2. Bounding Box Check ---
    float dot_product = (px - ax) * (bx - ax) + (py - ay) * (by - ay);
    
    return (dot_product >= -EPSILON_SQ_F) && (dot_product <= edge_len_sq + EPSILON_SQ_F);
}

// --- CUDA Kernel: Ray Casting (MODIFIED FOR MULTI-POLYGON) ---
__global__ void rayCastKernel(
    const float* value_x,  // [N]
    const float* value_y,  // [M]
    const float* xv,       // [V_clean_total] (Concatenated vertices)
    const float* yv,       // [V_clean_total]
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
    float px = value_x[c];
    float py = value_y[r];

    // --- 4. Initialize thread-local variables ---
    int totalCrossingCount = 0;
    bool totalOnEdge = false;

    // --- 5. Loop over all polygons ---
    for (size_t p = 0; p < numPolygons; ++p) {
        
        size_t poly_start_idx = poly_offsets[p];
        size_t poly_end_idx   = poly_offsets[p+1];
        size_t numVerticesInPoly = poly_end_idx - poly_start_idx;

        if (numVerticesInPoly < 4) { // e.g., v0, v1, v2, v0
            continue;
        }
        
        size_t numEdgesInPoly = numVerticesInPoly - 1; // It's pre-closed
        int polyCrossingCount = 0;
        bool polyOnEdge = false;

        // --- 5.A. Loop over all edges *in this polygon* ---
        for (size_t e_local = 0; e_local < numEdgesInPoly; ++e_local) {
            size_t e_global = poly_start_idx + e_local;
            
            float ax = xv[e_global];
            float ay = yv[e_global];
            float bx = xv[e_global + 1];
            float by = yv[e_global + 1];

            // --- 5.B. Check if the point is *on* this edge ---
            if (isPointOnSegment_device(px, py, ax, ay, bx, by)) {
                polyOnEdge = true;
                break; // No need to check other edges in *this* polygon
            }

            // --- 5.C. Ray-Casting Check ---
            bool y_span = ((ay <= py && py < by) || (by <= py && py < ay));
            
            if (y_span) {
                float x_intersect = (py - ay) * (bx - ax) / (by - ay) + ax;
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
            totalCrossingCount += polyCrossingCount;
        }
    } // End polygon loop

    // --- 6. Write final results ---
    onMesh[idx] = totalOnEdge;
    
    if (totalOnEdge) {
        inMesh[idx] = false; // 'on' takes precedence
    } else {
        inMesh[idx] = (totalCrossingCount % 2) != 0;
    }
}

// --- Main MEX Function (CPU Host) ---
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // --- 1. Check Inputs ---
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:nrhs", "Four inputs required: value_x, value_y, xv, yv.");
    }
    if (nlhs != 2) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:nlhs", "Two outputs required: [inMesh, onedgeMesh].");
    }
    
    // --- CHECK FOR SINGLE PRECISION ---
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]) || !mxIsSingle(prhs[2]) || !mxIsSingle(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:notSingle", "All inputs must be single precision (float).");
    }

    // --- 2. Get Grid Vectors (CPU) ---
    // --- FIX: Replaced mxGetSingles with mxGetData ---
    const float* value_x_cpu = (const float*)mxGetData(prhs[0]);
    const float* value_y_cpu = (const float*)mxGetData(prhs[1]);
    const size_t N = mxGetN(prhs[0]); // Columns (width)
    const size_t M = mxGetM(prhs[1]); // Rows (height)
    const size_t numGridPoints = M * N;

    // --- 3. Get Polygon Vectors (CPU) ---
    // --- FIX: Replaced mxGetSingles with mxGetData ---
    const float* xv_full_cpu = (const float*)mxGetData(prhs[2]);
    const float* yv_full_cpu = (const float*)mxGetData(prhs[3]);
    const size_t V_full = mxGetNumberOfElements(prhs[2]);
    if (V_full != mxGetNumberOfElements(prhs[3])) {
         mexErrMsgIdAndTxt("MATLAB:inpolygon_cuda:xyMismatch", "xv and yv must have the same number of elements.");
    }
    
    // --- 3.B. Pre-process Polygons (Handle NaNs) ---
    std::vector<size_t> nan_indices;
    std::vector<float> xv_clean_vec;
    std::vector<float> yv_clean_vec;
    std::vector<size_t> poly_offsets_vec;
    nan_indices.push_back(static_cast<size_t>(-1)); 
    for (size_t i = 0; i < V_full; ++i) {
        // Use std::isnan for floats
        if (std::isnan(xv_full_cpu[i]) || std::isnan(yv_full_cpu[i])) {
            nan_indices.push_back(i);
        }
    }
    nan_indices.push_back(V_full); 
    
    for (size_t p = 0; p < nan_indices.size() - 1; ++p) {
        size_t start_idx = nan_indices[p] + 1;
        size_t end_idx = nan_indices[p+1] - 1;
        
        if (end_idx < start_idx || (end_idx - start_idx + 1) < 2) continue;
        
        poly_offsets_vec.push_back(xv_clean_vec.size());
        
        for(size_t i = start_idx; i <= end_idx; ++i) {
            xv_clean_vec.push_back(xv_full_cpu[i]);
            yv_clean_vec.push_back(yv_full_cpu[i]);
        }
        
        // Check if the polygon is closed. If not, close it.
        float first_x = xv_clean_vec[poly_offsets_vec.back()];
        float first_y = yv_clean_vec[poly_offsets_vec.back()];
        float last_x  = xv_clean_vec.back();
        float last_y  = yv_clean_vec.back();
        if (std::abs(first_x - last_x) > EPSILON_F || std::abs(first_y - last_y) > EPSILON_F) {
            xv_clean_vec.push_back(first_x);
            yv_clean_vec.push_back(first_y);
        }
    }
    poly_offsets_vec.push_back(xv_clean_vec.size());
    
    const size_t numPolygons = poly_offsets_vec.size() - 1;
    const size_t V_clean_total = xv_clean_vec.size();
    if (numPolygons == 0 || V_clean_total == 0) {
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
    float* value_x_gpu;
    float* value_y_gpu;
    float* xv_gpu;
    float* yv_gpu;
    size_t* poly_offsets_gpu; 
    bool* inMesh_gpu;
    bool* onMesh_gpu;
    CUDA_CHECK(cudaMalloc(&value_x_gpu, N * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&value_y_gpu, M * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&xv_gpu, V_clean_total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&yv_gpu, V_clean_total * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&poly_offsets_gpu, (numPolygons + 1) * sizeof(size_t))); 
    CUDA_CHECK(cudaMalloc(&inMesh_gpu, numGridPoints * sizeof(bool)));
    CUDA_CHECK(cudaMalloc(&onMesh_gpu, numGridPoints * sizeof(bool)));

    // --- 6. Copy Data from CPU (Host) to GPU (Device) ---
    CUDA_CHECK(cudaMemcpy(value_x_gpu, value_x_cpu, N * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(value_y_gpu, value_y_cpu, M * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(xv_gpu, xv_clean_vec.data(), V_clean_total * sizeof(float), cudaMemcpyHostToDevice)); 
    CUDA_CHECK(cudaMemcpy(yv_gpu, yv_clean_vec.data(), V_clean_total * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(poly_offsets_gpu, poly_offsets_vec.data(), (numPolygons + 1) * sizeof(size_t), cudaMemcpyHostToDevice));
    
    CUDA_CHECK(cudaMemset(inMesh_gpu, 0, numGridPoints * sizeof(bool)));
    CUDA_CHECK(cudaMemset(onMesh_gpu, 0, numGridPoints * sizeof(bool)));

    // --- 7. Set up Kernel Launch Parameters ---
    int threadsPerBlock = 256; 
    int numBlocks = (numGridPoints + threadsPerBlock - 1) / threadsPerBlock;
    
    // --- 8. Launch the CUDA Kernel ---
    rayCastKernel<<<numBlocks, threadsPerBlock>>>(
        value_x_gpu, value_y_gpu,
        xv_gpu, yv_gpu,
        inMesh_gpu, onMesh_gpu,
        poly_offsets_gpu,    
        M, N,
        numPolygons           
    );
    
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    // --- 9. Copy Results from GPU (Device) to CPU (Host) ---
    CUDA_CHECK(cudaMemcpy(inMesh_cpu, inMesh_gpu, numGridPoints * sizeof(bool), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(onedgeMesh_cpu, onMesh_gpu, numGridPoints * sizeof(bool), cudaMemcpyDeviceToHost));

    // --- 10. Clean up GPU Memory ---
    CUDA_CHECK(cudaFree(value_x_gpu));
    CUDA_CHECK(cudaFree(value_y_gpu));
    CUDA_CHECK(cudaFree(xv_gpu));
    CUDA_CHECK(cudaFree(yv_gpu));
    CUDA_CHECK(cudaFree(poly_offsets_gpu));
    CUDA_CHECK(cudaFree(inMesh_gpu));
    CUDA_CHECK(cudaFree(onMesh_gpu));
}