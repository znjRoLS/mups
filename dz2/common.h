#define MAX 1024
#define ACCURACY 0.01

typedef struct {
    double time;
    double value;
} Result;

typedef struct {
    double time;
    int val_size;
    double *value;
} Result_Vect;

//int compare_and_print(double sequential_result, double parallel_result, double sequential_time, double parallel_time) {
//    int final_result = fabs(parallel_result - sequential_result) < ACCURACY;
//
//    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result - sequential_result);
//    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
//    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
//    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
//}

int compare_and_print(Result sequential_result, Result parallel_result, const char *label) {
    int final_result = fabs(parallel_result.value - sequential_result.value) < ACCURACY;

    printf("Results compare for %s:\n", label );
    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result.value - sequential_result.value);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_result.time - parallel_result.time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_result.time / parallel_result.time);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}

int compare_and_print_vect(Result_Vect sequential_result, Result_Vect parallel_result, const char *label) {

    int final_final_result = 1;

    printf("Results compare for %s:\n", label );

    if (sequential_result.val_size != parallel_result.val_size) {
        printf("Result size differs!! Sequental size %d, parallel size %d\n", sequential_result.val_size, parallel_result.val_size);
    }


    for (int i = 0 ; i < sequential_result.val_size; i ++) {
        int final_result = fabs(parallel_result.value[i] - sequential_result.value[i]) < ACCURACY;
        if (!final_result) {
            printf("Not accurate enough! sequental %f, parallel %f, diff %f, ind %d, %d\n",
                   sequential_result.value[i],
                   parallel_result.value[i],
                   fabs(parallel_result.value[i] - sequential_result.value[i]),
                   i/500, i % 500
            );
            final_final_result = 0;
        }
    }

    printf("Time gain (parallel) - (sequential) = %f\n", sequential_result.time - parallel_result.time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_result.time / parallel_result.time);
    printf( (final_final_result?"TEST PASSED\n":"TEST FAILED\n") );
}