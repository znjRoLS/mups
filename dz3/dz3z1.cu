# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

#include "common.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <device_functions.h>

//#define N 1024*1024*12

#define BLOCK_SIZE 256

#define MAX_BLOCKS 65535

#define MAX_MEMORY BLOCK_SIZE * MAX_BLOCKS

double f(double x) {
	double pi = 3.141592653589793;
	double value;

	value = 50.0 / (pi * (2500.0 * x * x + 1.0));

	return value;
}

__device__ double fDev(double x) {
	double pi = 3.141592653589793;
	double value;

	value = 50.0 / (pi * (2500.0 * x * x + 1.0));

	return value;
}

int sequential(int argc, char *argv[], Result_Vect *result) {
	double a;
	double b;
	double error;
	int i;
	int n;
	double total_q, total_t, total_s;
	double wtime_q, wtime_t, wtime_s;
	double x;
	double h;

	printf("\n\nSEQUENTIAL\n");
	
	result->time = 0;

	if (argc != 4) {
		n = 10000000;
		a = 0.0;
		b = 10.0;
	}
	else {
		n = atoi(argv[1]);
		a = atoi(argv[2]);
		b = atoi(argv[3]);
	}

	printf("\n");
	printf("QUAD:\n");
	printf("  Estimate the integral of f(x) from A to B.\n");
	printf("  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n");
	printf("\n");
	printf("  A        = %f\n", a);
	printf("  B        = %f\n", b);
	printf("  N        = %d\n", n);


	// Quadratic rule  
	wtime_q = omp_get_wtime();

	total_q = 0.0;

	for (i = 0; i < n; i++)
	{
		x = ((double)(n - i - 1) * a + (double)(i)* b) / (double)(n - 1);
		total_q = total_q + f(x);
	}

	wtime_q = omp_get_wtime() - wtime_q;

	total_q = (b - a) * total_q / (double)n;

	result->time += wtime_q;
	result->value[0] = total_q;


	// Trapezoidal rule  
	h = (b - a) / n;

	wtime_t = omp_get_wtime();

	total_t = 0.0;

	for (i = 0; i < n; i++)
	{
		x = a + i * h;
		if (i > 0 && i < n - 1)
			total_t = total_t + f(x);
		else
			total_t = total_t + 0.5 * f(x);
	}

	total_t = h * total_t;

	wtime_t = omp_get_wtime() - wtime_t;


	result->time += wtime_t;
	result->value[1] = total_t;

	// Simpson 1/3 rule  

	h = (b - a) / n;

	wtime_s = omp_get_wtime();

	total_s = 0.0;

	for (i = 0; i < n; i++)
	{
		x = a + i * h;
		if (i == 0 || i == n - 1)
			total_s = total_s + f(x);
		else if (i % 2 == 1)
			total_s = total_s + 4 * f(x);
		else
			total_s = total_s + 2 * f(x);
	}

	total_s = h / 3 * total_s;

	wtime_s = omp_get_wtime() - wtime_s;

	result->time += wtime_s;
	result->value[2] = total_s;

	printf("\n");
	printf("  Estimate quadratic rule = %24.16f\n", total_q);
	printf("  Estimate trapezoidal rule = %24.16f\n", total_t);
	printf("  Estimate Simpson 1/3 rule = %24.16f\n", total_s);
	printf("  Time quadratic rule = %f\n", wtime_q);
	printf("  Time trapezoidal rule = %f\n", wtime_t);
	printf("  Time Simpson 1/3 rule = %f\n", wtime_s);
	printf("\n");
	printf("  Normal end of execution.\n");
	printf("\n");

	return 0;
}

/////////////////// parallel

// Simple reduction kernel
__global__ void reductionSumKernel(double* devA, double* blockResults, int n) {
	extern __shared__ double sharedData[];

	unsigned int tid = threadIdx.x;
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;

	// Load block in the shared memory
	if (i < n) sharedData[tid] = devA[i];
	else sharedData[tid] = 0;

	__syncthreads();

	// Do reduction in shared memory
	for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
		if (tid < s) {
			sharedData[tid] += sharedData[tid + s];
		}
		__syncthreads();
	}

	// Write result for this block to global memory 
	if (tid == 0) blockResults[blockIdx.x] = sharedData[0];
}

double sumReduction(double* devA, int n) {
	double gpuSum = 0;
	int numBlocks = 0;
	double *devBlockRes;
	
	// Run kernel several times until the work is done
	numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
	//printf("this is numblocks %d\n", numBlocks);
	cudaMalloc((void **)&devBlockRes, numBlocks * sizeof(double));

	reductionSumKernel <<< numBlocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double) >>>(devA, devBlockRes, n);

	while (numBlocks > 1) {
		n = numBlocks;
		numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
		reductionSumKernel <<< numBlocks, BLOCK_SIZE, BLOCK_SIZE * sizeof(double) >>>(devBlockRes, devBlockRes, n);
	}

	// Copy back the results
	cudaMemcpy(&gpuSum, devBlockRes, sizeof(double), cudaMemcpyDeviceToHost);

	cudaFree(devBlockRes);

	return gpuSum;
}

__global__ void compute_kernel_quad_big(double *devA, double a, double b, int n, int offset, int nSize) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x + offset;
	while (i < nSize + offset) {
		double x = ((double)(n - i - 1) * a + (double)(i)* b) / (double)(n - 1);
		devA[i-offset] = fDev(x);
		i += blockDim.x * gridDim.x;
	}
}


__global__ void compute_kernel_trapezoidal_big(double *devA, double a, double b, int n, int offset, int nSize) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x + offset;
	while (i < nSize + offset) {
		double h = (b - a) / n;

		double x = a + i * h;
		if (i > 0 && i < n - 1)
			devA[i - offset] = fDev(x);
		else
			devA[i - offset] = 0.5 * fDev(x);

		i += blockDim.x * gridDim.x;
	}
}

__global__ void compute_kernel_simpson_big(double *devA, double a, double b, int n, int offset, int nSize) {
	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x + offset;
	if (i < nSize + offset) {
		double h = (b - a) / n;

		double x = a + i * h;
		if (i == 0 || i == n - 1)
			devA[i-offset] = fDev(x);
		else if (i % 2 == 1)
			devA[i-offset] = 4 * fDev(x);
		else
			devA[i-offset] = 2 * fDev(x);

		i += blockDim.x * gridDim.x;
	}
}

//__global__ void compute_kernel_quad(double *devA, double a, double b, int n) {
//	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i < n) {
//		double x = ((double)(n - i - 1) * a + (double)(i)* b) / (double)(n - 1);
//		devA[i] = fDev(x);
//	}
//}

//__global__ void compute_kernel_trapezoidal(double *devA, double a, double b, int n) {
//	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i < n) {
//		double h = (b - a) / n;
//
//		double x = a + i * h;
//		if (i > 0 && i < n - 1)
//			devA[i] = fDev(x);
//		else
//			devA[i] = 0.5 * fDev(x);
//	}
//}

//double parallel_trapezoidal(double a, double b, int n) {
//	double *devA;
//
//	cudaMalloc((void **)&devA, n * sizeof(double));
//	//cudaMemcpy(devA, A, size, cudaMemcpyHostToDevice);
//
//	int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
//
//	compute_kernel_trapezoidal << < numBlocks, BLOCK_SIZE >> > (devA, a, b, n);
//	double total_q = sumReduction(devA, n);
//
//	cudaFree(devA);
//
//	return total_q;
//}

//__global__ void compute_kernel_simpson(double *devA, double a, double b, int n) {
//	unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if (i < n) {
//		double h = (b - a) / n;
//
//		double x = a + i * h;		
//		if (i == 0 || i == n - 1)
//			devA[i] = fDev(x);
//		else if (i % 2 == 1)
//			devA[i] = 4 * fDev(x);
//		else
//			devA[i] = 2 * fDev(x);
//	}
//}

//double parallel_simpson(double a, double b, int n) {
//	double *devA;
//
//	cudaMalloc((void **)&devA, n * sizeof(double));
//	//cudaMemcpy(devA, A, size, cudaMemcpyHostToDevice);
//
//	int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
//
//	compute_kernel_simpson << < numBlocks, BLOCK_SIZE >> > (devA, a, b, n);
//	double total_q = sumReduction(devA, n);
//
//	cudaFree(devA);
//
//	return total_q;
//}


double parallel_compute(double a, double b, int n, void (*kernel)(double*,double,double,int,int,int)) {
	double *devA;
	double total_q = 0;

	for (int ni = 0; ni < n; ni += MAX_MEMORY) {
		int nSize = MAX_MEMORY;
		if (n - ni < MAX_MEMORY)
			nSize = n - ni;

		cudaMalloc((void **)&devA, nSize * sizeof(double));
		if (cudaSuccess != cudaGetLastError()) {
			printf("couldnt allocate %d doubles\n", nSize);
			break;
		}

		int numBlocks = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
		if (numBlocks > MAX_BLOCKS)
			numBlocks = MAX_BLOCKS;

		kernel << < numBlocks, BLOCK_SIZE >> > (devA, a, b, n, ni, nSize);
		if (cudaSuccess != cudaGetLastError()) {
			printf("error computing\n");
			break;
		}

		total_q += sumReduction(devA, nSize);
		if (cudaSuccess != cudaGetLastError()) {
			printf("error reducting\n");
			break;
		}
		
		cudaFree(devA);
	}

	return total_q;	
}

int parallel(int argc, char *argv[], Result_Vect *result) {
	double a;
	double b;
	double error;
	int i;
	int n;
	double total_q, total_t, total_s;
	double wtime_q, wtime_t, wtime_s;
	double x;
	double h;

	printf("\n\nPARALLEL\n");
	
	result->time = 0;

	if (argc != 4) {
		n = 10000000;
		a = 0.0;
		b = 10.0;
	}
	else {
		n = atoi(argv[1]);
		a = atoi(argv[2]);
		b = atoi(argv[3]);
	}

	printf("\n");
	printf("QUAD:\n");
	printf("  Estimate the integral of f(x) from A to B.\n");
	printf("  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n");
	printf("\n");
	printf("  A        = %f\n", a);
	printf("  B        = %f\n", b);
	printf("  N        = %d\n", n);


	// Quadratic rule  
	wtime_q = omp_get_wtime();

	/*total_q = 0.0;

	for (i = 0; i < n; i++)
	{
		x = ((double)(n - i - 1) * a + (double)(i)* b) / (double)(n - 1);
		total_q = total_q + f(x);
	}*/

	total_q = parallel_compute(a, b, n, compute_kernel_quad_big);

	total_q = (b - a) * total_q / (double)n;

	wtime_q = omp_get_wtime() - wtime_q;



	result->time += wtime_q;
	result->value[0] = total_q;

	// Trapezoidal rule  
	h = (b - a) / n;

	wtime_t = omp_get_wtime();

	/*total_t = 0.0;

	for (i = 0; i < n; i++)
	{
		x = a + i * h;
		if (i > 0 && i < n - 1)
			total_t = total_t + f(x);
		else
			total_t = total_t + 0.5 * f(x);
	}*/

	/*total_t = parallel_trapezoidal(a, b, n);*/
	total_t = parallel_compute(a, b, n, compute_kernel_trapezoidal_big);

	total_t = h * total_t;

	wtime_t = omp_get_wtime() - wtime_t;

	result->time += wtime_t;
	result->value[1] = total_t;

	// Simpson 1/3 rule  

	h = (b - a) / n;

	wtime_s = omp_get_wtime();
/*
	total_s = 0.0;

	for (i = 0; i < n; i++)
	{
		x = a + i * h;
		if (i == 0 || i == n - 1)
			total_s = total_s + f(x);
		else if (i % 2 == 1)
			total_s = total_s + 4 * f(x);
		else
			total_s = total_s + 2 * f(x);
	}*/

	total_s = parallel_compute(a, b, n, compute_kernel_simpson_big);
	/*total_s = parallel_simpson(a, b, n);*/

	total_s = h / 3 * total_s;

	wtime_s = omp_get_wtime() - wtime_s;

	result->time += wtime_s;
	result->value[2] = total_s;

	printf("\n");
	printf("  Estimate quadratic rule = %24.16f\n", total_q);
	printf("  Estimate trapezoidal rule = %24.16f\n", total_t);
	printf("  Estimate Simpson 1/3 rule = %24.16f\n", total_s);
	printf("  Time quadratic rule = %f\n", wtime_q);
	printf("  Time trapezoidal rule = %f\n", wtime_t);
	printf("  Time Simpson 1/3 rule = %f\n", wtime_s);
	printf("\n");
	printf("  Normal end of execution.\n");
	printf("\n");

	return 0;
}


int main(int argc, char *argv[]) {

	//double sequential_result, parallel_result, sequential_time, parallel_time;

	Result_Vect seq_result;
	Result_Vect par_result;

	seq_result.val_size = 3;
	seq_result.value = (double*)malloc(3 * sizeof(double));
	par_result.val_size = 3;
	par_result.value = (double*)malloc(3 * sizeof(double));


	/*for (int i = 1; ; i <<= 1) {
		double *nekid;
		cudaMalloc((void **)&nekid, i * sizeof(double));
		if (cudaSuccess != cudaGetLastError()) {
			printf("couldnt allocate %d doubles\n", i);
			break;
		}
		else {
			printf("allocated %d doubles\n", i);
			cudaFree(nekid);
		}

	}*/

	sequential(argc, argv, &seq_result);
	parallel(argc, argv, &par_result);

	compare_and_print_vect(seq_result, par_result, "Numeric integration");
}
