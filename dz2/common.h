#define MAX 1024
#define ACCURACY 0.01

typedef struct {
    double time;
    double value;
} Result;

int compare_and_print(double sequential_result, double parallel_result, double sequential_time, double parallel_time) {
    int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result - sequential_result);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}

int compare_and_print_result(Result sequential_result, Result parallel_result, const char *label) {
    int final_result = fabs(parallel_result.value - sequential_result.value) < ACCURACY;

    printf("Results compare for %s:\n", label );
    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result.value - sequential_result.value);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_result.time - parallel_result.time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_result.time / parallel_result.time);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}