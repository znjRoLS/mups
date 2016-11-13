#define MAX 1024
#define ACCURACY 0.01
#define MAXTHREADS 8

#include <math.h>
#include <stdio.h>

int finish_1(double sequential_result, double parallel_result, double sequential_time, double parallel_time) {
    int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    printf("Result difference (parallel) - (sequential) = %f\n", parallel_result - sequential_result);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
    printf( (final_result?"TEST PASSED\n":"TEST FAILED\n") );
}

int finish_2(
        unsigned char * sequential_result,
        unsigned char * parallel_result,
        unsigned int sequential_size,
        unsigned int parallel_size,
        double sequential_time,
        double parallel_time
) {
    //int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    int final_result = 0;
    if (sequential_size != parallel_size) {
        printf("Sizes aren't equal!\n");
        printf("TEST FAILED\n");
    }

    for (int i = 0 ; i < sequential_size; i ++) {
        final_result += abs(sequential_result[i] - parallel_result[i]);
        if (abs(sequential_result[i] - parallel_result[i]) != 0) {
            printf ("diff at %d\n", i);
        }
    }

    printf("\n----------------------------------\n");
    printf("Result difference (parallel) - (sequential) = %d\n", final_result);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
    printf( (final_result?"TEST FAILED\n":"TEST PASSED\n") );
    printf("\n----------------------------------\n");
}

int finish_3(
        int * sequential_result,
        int * parallel_result,
        unsigned int sequential_size,
        unsigned int parallel_size,
        double sequential_time,
        double parallel_time
) {
    //int final_result = fabs(parallel_result - sequential_result) < ACCURACY;

    int final_result = 0;
    if (sequential_size != parallel_size) {
        printf("Sizes aren't equal!\n");
        printf("TEST FAILED\n");
    }

    for (int i = 0 ; i < sequential_size; i ++) {
        final_result += abs(sequential_result[i] - parallel_result[i]);
        if (abs(sequential_result[i] - parallel_result[i]) != 0) {
            printf ("diff at %d %d %d\n", i, sequential_result[i], parallel_result[i]);
        }
    }

    printf("\n----------------------------------\n");
    printf("Result difference (parallel) - (sequential) = %d\n", final_result);
    printf("Time gain (parallel) - (sequential) = %f\n", sequential_time - parallel_time);
    printf("Time acceleration (sequential)/(parallel) = %f\n", sequential_time / parallel_time);
    printf( (final_result?"TEST FAILED\n":"TEST PASSED\n") );
    printf("\n----------------------------------\n");
}