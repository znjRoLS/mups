# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

#include "common.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <device_launch_parameters.h>
#include <device_functions.h>


double cpu_time(void)
{
	double value;

	value = (double)clock() / (double)CLOCKS_PER_SEC;

	return value;
}


int sequential(int argc, char *argv[], Result_Vect *result)
{
	int M;
	int N;

	double ctime;
	double ctime1;
	double ctime2;
	double diff;
	double epsilon;
	FILE *fp;
	int i;
	int iterations;
	int iterations_print;
	int j;
	double mean;
	char output_file[80];
	int success;

	double **u;
	double **w;

	
	printf("\n\nSEQUENTIAL\n");
	
	if (argc != 5) {
		printf("Wrong number of arguments!\n");
		return 1;
	}
	else {
		success = sscanf(argv[1], "%d", &M);
		success += sscanf(argv[2], "%d", &N);
		success += sscanf(argv[3], "%lf", &epsilon);
		success += sscanf(argv[4], "%s", output_file);

		if (success != 4) {
			printf("Wrong arguments!\n");
			return 2;
		}
	}

	printf("\n");
	printf("HEATED_PLATE\n");
	printf("  C version\n");
	printf("  A program to solve for the steady state temperature distribution\n");
	printf("  over a rectangular plate.\n");
	printf("\n");
	printf("  Spatial grid of %d by %d points.\n", M, N);
	printf("\n");
	printf("  The iteration will be repeated until the change is <= %f\n", epsilon);
	diff = epsilon;
	printf("\n");
	printf("  The steady state solution will be written to %s.\n", output_file);

	u = (double **)malloc(M * sizeof(double*));
	for (i = 0; i < M; i++)
		u[i] = (double *)malloc(N * sizeof(double));

	w = (double **)malloc(M * sizeof(double*));
	for (i = 0; i < M; i++)
		w[i] = (double *)malloc(N * sizeof(double));

	/*
	Set the boundary values, which don't change.
	*/
	for (i = 1; i < M - 1; i++)
	{
		w[i][0] = 100.0;
	}
	for (i = 1; i < M - 1; i++)
	{
		w[i][N - 1] = 100.0;
	}
	for (j = 0; j < N; j++)
	{
		w[M - 1][j] = 100.0;
	}
	for (j = 0; j < N; j++)
	{
		w[0][j] = 0.0;
	}
	/*
	Average the boundary values, to come up with a reasonable
	initial value for the interior.
	*/
	mean = 0.0;
	for (i = 1; i < M - 1; i++)
	{
		mean = mean + w[i][0];
	}
	for (i = 1; i < M - 1; i++)
	{
		mean = mean + w[i][N - 1];
	}
	for (j = 0; j < N; j++)
	{
		mean = mean + w[M - 1][j];
	}
	for (j = 0; j < N; j++)
	{
		mean = mean + w[0][j];
	}
	mean = mean / (double)(2 * M + 2 * N - 4);
	/*
	Initialize the interior solution to the mean value.
	*/
	for (i = 1; i < M - 1; i++)
	{
		for (j = 1; j < N - 1; j++)
		{
			w[i][j] = mean;
		}
	}
	/*
	iterate until the  new solution W differs from the old solution U
	by no more than EPSILON.
	*/
	iterations = 0;
	iterations_print = 1;
	printf("\n");
	printf(" Iteration  Change\n");
	printf("\n");
	ctime1 = cpu_time();

	while (epsilon <= diff)
	{
		/*
		Save the old solution in U.
		*/
		for (i = 0; i < M; i++)
		{
			for (j = 0; j < N; j++)
			{
				u[i][j] = w[i][j];
			}
		}
		/*
		Determine the new estimate of the solution at the interior points.
		The new solution W is the average of north, south, east and west neighbors.
		*/
		diff = 0.0;


		for (i = 1; i < M - 1; i++)
		{
			for (j = 1; j < N - 1; j++)
			{
				w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;

				if (diff < fabs(w[i][j] - u[i][j]))
				{
					diff = fabs(w[i][j] - u[i][j]);
				}
			}
		}

		/*for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
				if (i < 2 && j < 2) {
					printf("seq - devvv? %d %d %llf - %llf\n", i, j, u[i][j], w[i][j]);
				}*/

		iterations++;
		if (iterations == iterations_print)
		{
			printf("  %8d  %f\n", iterations, diff);
			iterations_print = 2 * iterations_print;
		}

		/*if (iterations == 5) {
			diff = 0;
		}*/
	}
	ctime2 = cpu_time();
	ctime = ctime2 - ctime1;

	printf("\n");
	printf("  %8d  %f\n", iterations, diff);
	printf("\n");
	printf("  Error tolerance achieved.\n");
	printf("  CPU time = %f\n", ctime);
	/*
	Write the solution to the output file.
	*/
	fp = fopen(output_file, "w");

	fprintf(fp, "%d\n", M);
	fprintf(fp, "%d\n", N);

	result->val_size = M*N;
	result->value = (double*)malloc(M*N * sizeof(double));
	result->time = ctime;

	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			fprintf(fp, "%6.2f ", w[i][j]);
			result->value[i*N + j] = w[i][j];
		}
		fputc('\n', fp);
	}
	fclose(fp);

	printf("\n");
	printf("  Solution written to the output file %s\n", output_file);
	/*
	All done!
	*/
	printf("\n");
	printf("HEATED_PLATE:\n");
	printf("  Normal end of execution.\n");

	return 0;

}

////////////////////////// parallel

//__device__ static double atomicMax(double* address, double val)
//{
//	unsigned long long int* address_as_i = 
//		(unsigned long long int*)address;
//	unsigned long long int old = *address_as_i, assumed;
//	do {
//		assumed = old;
//		old = ::atomicCAS(address_as_i, assumed,
//			__double_as_longlong( (val > __longlong_as_double(assumed) ) ? val : __longlong_as_double(assumed)));
//	} while (assumed != old);
//	return __longlong_as_double(old);
//}

__device__ void atomicMax(double * const address, const double value)
{
	if (*address >= value)
	{
		return;
	}

	unsigned long long int * const address_as_i = (unsigned long long int *)address;
	unsigned long long int old = *address_as_i, assumed;

	do
	{
		assumed = old;
		if (__longlong_as_double(assumed) >= value)
		{
			break;
		}

		old = atomicCAS(address_as_i, assumed, __double_as_longlong(value));
	} while (assumed != old);
}

__global__ void heated_kernel(double *devA, double *devB, int N, int M, double *epsilon) {
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;


	//printf("nesto? %d %d\n", i ,j);
	if (i > 0 && j > 0 && i < M-1 && j < N-1) {
		devB[i*N + j] = (devA[(i - 1)*N + j] + devA[(i + 1)*N + j] + devA[i * N + j - 1] + devA[i*N + j + 1]) / 4.0;
		
		atomicMax(epsilon, devB[i*N + j] - devA[i*N + j]);
		atomicMax(epsilon, devA[i*N + j] - devB[i*N + j]);
	}

	/*if (i < 2 && j < 2) {
		printf("devvv? %d %d %llf - %llf\n", i, j, devA[i*N + j], devB[i*N + j]);
	}*/
}

__global__ void heated_kernel2(double *devA, double *devB, int N, int M, double epsilon, int *isBigger) {
	int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
	int j = blockIdx.y * blockDim.y + threadIdx.y + 1;


	//printf("nesto? %d %d\n", i ,j);
	if (i > 0 && j > 0 && i < M - 1 && j < N - 1) {
		devB[i*N + j] = (devA[(i - 1)*N + j] + devA[(i + 1)*N + j] + devA[i * N + j - 1] + devA[i*N + j + 1]) / 4.0;

		if (devB[i*N + j] - devA[i*N + j] > epsilon) {
			*isBigger = 1;
		}
		if (devA[i*N + j] - devB[i*N + j] > epsilon) {
			*isBigger = 1;
		}
	}

	/*if (i < 2 && j < 2) {
	printf("devvv? %d %d %llf - %llf\n", i, j, devA[i*N + j], devB[i*N + j]);
	}*/
}

#define USEPLATE 1

void parallel_heated_plate(double **u, double **w, int N, int M, double epsilon) {

	double diff = epsilon;
	int iterations = 0;
	int iterations_print = 1;
	int i, j;
	double zeroDouble = 0.0;
	int zeroInt = 0;
	int isDiff = 1;

	dim3 threadsPerBlock(32, 16);
	dim3 numBlocks((M - 2 +threadsPerBlock.x - 1) / threadsPerBlock.x, (N - 2 + threadsPerBlock.y - 1) / threadsPerBlock.y);

	double *devU, *devW;
	double *epsilonCuda;
	int *isDiffCuda;

	cudaMalloc((void **)&devU, N*M * sizeof(double));
	cudaMalloc((void **)&devW, N*M * sizeof(double));
	if (cudaSuccess != cudaGetLastError()) {
		printf("errorr------------------------------------------\n");
	}
	cudaMalloc((void **)&epsilonCuda, sizeof(double));
	cudaMalloc((void **)&isDiffCuda, sizeof(int));

	if (cudaSuccess != cudaGetLastError()) {
		printf("errorr------------------------------------------\n");
	}

	for (int i = 0; i < M; i++) {
		cudaMemcpy(devU + i*N, u[i], N * sizeof(double), cudaMemcpyHostToDevice);

		if (cudaSuccess != cudaGetLastError()) {
			printf("errorr------------------------------------------\n");
		}
	}

	for (int i = 0; i < M; i++) {
		cudaMemcpy(devW + i*N, w[i], N * sizeof(double), cudaMemcpyHostToDevice);

		if (cudaSuccess != cudaGetLastError()) {
			printf("errorr------------------------------------------\n");
		}
	}

	//cudaMemcpy(devU, u, N*M * sizeof(double), cudaMemcpyHostToDevice);

	int cnt = 0;

#if USEPLATE == 1
	while (epsilon <= diff)
#endif	
#if USEPLATE == 2
	while (isDiff)
#endif	
	{
		cnt++;

		if (USEPLATE == 1) {
			cudaMemcpy(epsilonCuda, &zeroDouble, sizeof(double), cudaMemcpyHostToDevice);
		}
		else if (USEPLATE == 2) {
			cudaMemcpy(isDiffCuda, &zeroInt, sizeof(int), cudaMemcpyHostToDevice);
		}
		
		/*
		Determine the new estimate of the solution at the interior points.
		The new solution W is the average of north, south, east and west neighbors.
		*/
		if (cnt % 2) {
			/*printf("%d u w\n", cnt);*/
			if (USEPLATE == 1) {
				heated_kernel << < numBlocks, threadsPerBlock >> > (devU, devW, N, M, epsilonCuda);
			}
			else if (USEPLATE == 2) {
				heated_kernel2 << < numBlocks, threadsPerBlock >> > (devU, devW, N, M, epsilon, isDiffCuda);
			}
			
		}
		else {
			/*printf("%d w u\n", cnt);*/
			
			if (USEPLATE == 1) {
				heated_kernel << < numBlocks, threadsPerBlock >> > (devW, devU, N, M, epsilonCuda);
			}
			else if (USEPLATE == 2) {
				heated_kernel2 << < numBlocks, threadsPerBlock >> > (devW, devU, N, M, epsilon, isDiffCuda);
			}
			
		}

		cudaDeviceSynchronize();
		if (cudaSuccess != cudaGetLastError()) {
			printf("errorr------------------------------------------\n");
		}

		if (USEPLATE == 1) {
			cudaMemcpy(&diff, epsilonCuda, sizeof(double), cudaMemcpyDeviceToHost);
		}
		else if (USEPLATE == 2) {
			cudaMemcpy(&isDiff, isDiffCuda, sizeof(int), cudaMemcpyDeviceToHost);
		}
		
		//if (cudaSuccess != cudaGetLastError()) {
		//	printf("errorr\n");
		//}

		iterations++;
		if (iterations == iterations_print)
		{
			printf("  %8d  %f\n", iterations, diff);
			iterations_print = 2 * iterations_print;
		}
		/*if (iterations == 5) {
			diff = 0;
		}*/
	}

	if (cnt % 2) {
		for (int i = 0; i < M; i++) {
			cudaMemcpy(w[i], devW + i*N, N * sizeof(double), cudaMemcpyDeviceToHost);
			if (cudaSuccess != cudaGetLastError()) {
				printf("errorr\n");
			}
			//printf("sta bre? %d, %llf, %llf, %llf\n", i, w[i][0], w[i][1], w[i][2]);
		}
	}
	else {
		for (int i = 0; i < M; i++) {
			cudaMemcpy(w[i], devU + i*N, N * sizeof(double), cudaMemcpyDeviceToHost);
			if (cudaSuccess != cudaGetLastError()) {
				printf("errorr\n");
			}
			//printf("sta bru? %d, %llf, %llf, %llf\n", i, w[i][0], w[i][1], w[i][2]);
		}
	}

	printf("\n");
	printf("  %8d  %f\n", iterations, diff);

}


int parallel(int argc, char *argv[], Result_Vect *result)
{
	int M;
	int N;

	double ctime;
	double ctime1;
	double ctime2;
	double diff;
	double epsilon;
	FILE *fp;
	int i;
	int iterations;
	int iterations_print;
	int j;
	double mean;
	char output_file[80];
	int success;

	double **u;
	double **w;

	
	printf("\n\nPARALLEL\n");
	
	if (argc != 5) {
		printf("Wrong number of arguments!\n");
		return 1;
	}
	else {
		success = sscanf(argv[1], "%d", &M);
		success += sscanf(argv[2], "%d", &N);
		success += sscanf(argv[3], "%lf", &epsilon);
		success += sscanf(argv[4], "%s", output_file);

		if (success != 4) {
			printf("Wrong arguments!\n");
			return 2;
		}
	}

	printf("\n");
	printf("HEATED_PLATE\n");
	printf("  C version\n");
	printf("  A program to solve for the steady state temperature distribution\n");
	printf("  over a rectangular plate.\n");
	printf("\n");
	printf("  Spatial grid of %d by %d points.\n", M, N);
	printf("\n");
	printf("  The iteration will be repeated until the change is <= %f\n", epsilon);
	diff = epsilon;
	printf("\n");
	printf("  The steady state solution will be written to %s.\n", output_file);

	u = (double **)malloc(M * sizeof(double*));
	for (i = 0; i < M; i++)
		u[i] = (double *)malloc(N * sizeof(double));

	w = (double **)malloc(M * sizeof(double*));
	for (i = 0; i < M; i++)
		w[i] = (double *)malloc(N * sizeof(double));

	/*
	Set the boundary values, which don't change.
	*/
	for (i = 1; i < M - 1; i++)
	{
		w[i][0] = 100.0;
	}
	for (i = 1; i < M - 1; i++)
	{
		w[i][N - 1] = 100.0;
	}
	for (j = 0; j < N; j++)
	{
		w[M - 1][j] = 100.0;
	}
	for (j = 0; j < N; j++)
	{
		w[0][j] = 0.0;
	}
	/*
	Average the boundary values, to come up with a reasonable
	initial value for the interior.
	*/
	mean = 0.0;
	for (i = 1; i < M - 1; i++)
	{
		mean = mean + w[i][0];
	}
	for (i = 1; i < M - 1; i++)
	{
		mean = mean + w[i][N - 1];
	}
	for (j = 0; j < N; j++)
	{
		mean = mean + w[M - 1][j];
	}
	for (j = 0; j < N; j++)
	{
		mean = mean + w[0][j];
	}
	mean = mean / (double)(2 * M + 2 * N - 4);
	/*
	Initialize the interior solution to the mean value.
	*/
	for (i = 1; i < M - 1; i++)
	{
		for (j = 1; j < N - 1; j++)
		{
			w[i][j] = mean;
		}
	}
	for (i = 0; i < M; i++)
	{
		for (j = 0; j < N; j++)
		{
			u[i][j] = w[i][j];
		}
	}

	/*
	iterate until the  new solution W differs from the old solution U
	by no more than EPSILON.
	*/
	iterations = 0;
	iterations_print = 1;
	printf("\n");
	printf(" Iteration  Change\n");
	printf("\n");
	ctime1 = cpu_time();

	parallel_heated_plate(u, w, N, M, epsilon);

	//while (epsilon <= diff)
	//{
	//	/*
	//	Save the old solution in U.
	//	*/
	//	for (i = 0; i < M; i++)
	//	{
	//		for (j = 0; j < N; j++)
	//		{
	//			u[i][j] = w[i][j];
	//		}
	//	}
	//	/*
	//	Determine the new estimate of the solution at the interior points.
	//	The new solution W is the average of north, south, east and west neighbors.
	//	*/
	//	diff = 0.0;
	//	for (i = 1; i < M - 1; i++)
	//	{
	//		for (j = 1; j < N - 1; j++)
	//		{
	//			w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;

	//			if (diff < fabs(w[i][j] - u[i][j]))
	//			{
	//				diff = fabs(w[i][j] - u[i][j]);
	//			}
	//		}
	//	}
	//	iterations++;
	//	if (iterations == iterations_print)
	//	{
	//		printf("  %8d  %f\n", iterations, diff);
	//		iterations_print = 2 * iterations_print;
	//	}
	//}
	ctime2 = cpu_time();
	ctime = ctime2 - ctime1;

	printf("\n");
	printf("  Error tolerance achieved.\n");
	printf("  CPU time = %f\n", ctime);
	/*
	Write the solution to the output file.
	*/
	fp = fopen(output_file, "w");

	fprintf(fp, "%d\n", M);
	fprintf(fp, "%d\n", N);

	result->val_size = M*N;
	result->value = (double*)malloc(M*N * sizeof(double));
	result->time = ctime;

	for (i = 0; i < M; i++)
	{
		//printf("sta brej? %d, %llf, %llf, %llf\n", i, w[i][0], w[i][1], w[i][2]);

		for (j = 0; j < N; j++)
		{
			fprintf(fp, "%6.2f ", w[i][j]);
			result->value[i*N + j] = w[i][j];
		}
		fputc('\n', fp);
	}
	fclose(fp);

	printf("\n");
	printf("  Solution written to the output file %s\n", output_file);
	/*
	All done!
	*/
	printf("\n");
	printf("HEATED_PLATE:\n");
	printf("  Normal end of execution.\n");

	return 0;

}

int main(int argc, char * argv[]) {
	Result_Vect seq_result, par_result;

	sequential(argc, argv, &seq_result);
	parallel(argc, argv, &par_result);

	compare_and_print_vect(seq_result, par_result, "heated plate");
}