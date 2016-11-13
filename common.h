#define MAX 1024
#define ACCURACY 0.01

int finish(double sequential_result, double parallel_result) {
    int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    printf("Difference (parallel) - (sequential) = %f\n", parallel_result - sequential_result);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}