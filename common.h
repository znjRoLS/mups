#define MAX 1024
#define ACCURACY 0.01

int finish(double sequential_result, double parallel_result, double sequential_time, double parallel_time) {
    int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result - sequential_result);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}