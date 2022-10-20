
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#include "../deps/CudaCommon/helper_math.h"

cudaError_t P2G_with_Cuda(// PARTICLES
    double2* particles_x, double pm,

    // GRID
    double* grid_m, double grid_dx, double2 grid_min_point, int grid_max_i, int grid_max_j);


__device__ double CubicBSpline(double x)
{
    x = abs(x);
    if (0.0 <= x && x < 1.0) {
        return 0.5 * x * x * x - x * x + 2.0 / 3.0;
    }
    else if (1.0 <= x && x < 2.0) {
        return (2.0 - x) * (2.0 - x) * (2.0 - x) / 6.0;
    }
    else {
        return 0.0;
    }
}

__global__ void P2G_Kernel(
    // PARTICLES
    double2* particles_x, double pm, 
    
    // GRID
    double* grid_m, double grid_dx, double2 grid_min_point, int grid_max_i, int grid_max_j)
{
    // Assume grid_m has been set to 0 for all indices


    int p_id = threadIdx.x;

    double2 px = particles_x[p_id];

    double2 relative_point;
    relative_point.x = px.x - grid_min_point.x;
    relative_point.y = px.y - grid_min_point.y;

    int bot_left_index_i = (int)floor(relative_point.x / grid_dx) - 1;
    int bot_left_index_j = (int)floor(relative_point.y / grid_dx) - 1;

    // CUBIC B-SPLINE
    for (int i = 0; i <= 3; i++) 
    {
        for (int j = 0; j <= 3; j++) 
        {
            // conditional branching... thats ok right?
            bool in_bounds = 0 < bot_left_index_i && bot_left_index_i < grid_max_i &&
                0 < bot_left_index_j && bot_left_index_j < grid_max_j;
            

            if (in_bounds)
            {
                grid_m[bot_left_index_i * grid_max_j + bot_left_index_j] += pm;
            }
        }
    }
}

__global__ void Grid_Reset_Kernel(double* grid_m, int grid_max_i, int grid_max_j)
{
    int i = threadIdx.x;
    int j = threadIdx.y;
    grid_m[i * grid_max_j + j] = 0;
}



int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to project particle masses onto the grid
cudaError_t P2G_with_Cuda(
    // PARTICLES
    double2* particles_x, double pm,

    // GRID
    double* grid_m, double grid_dx, double2 grid_min_point, int grid_max_i, int grid_max_j)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
