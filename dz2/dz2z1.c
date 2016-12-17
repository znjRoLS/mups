//
// Created by rols on 12/17/16.
//

# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

#include <mpi.h>
#include <omp.h>
#include "common.h"

int myrank, commSize;

double f ( double x );

double process(double a, double b, int n) {
    double total = 0.0;
    int i;
    double x;

    for ( i = 0; i < n; i++ )
    {
        x = ( ( double ) ( n - i - 1 ) * a + ( double ) ( i ) * b ) / ( double ) ( n - 1 );
        total = total + f ( x );
    }

    total = ( b - a ) * total / ( double ) n;

    return total;
}

int parallel ( int argc, char *argv[] , double *result, double *time) {
    double a;
    double b;
    double error;
    double exact = 0.49936338107645674464;
    int i;
    int n;
    double total;
    double wtime;

    double a_my, b_my, n_my, chunkSize, total_my;

    printf("PARALLEL RUN\n");


    if (argc != 4) {
        n = 10000000;
        a = 0.0;
        b = 10.0;
    } else {
        n = atoi(argv[1]);
        a = atoi(argv[2]);
        b = atoi(argv[3]);
    }

    if (myrank == 0) {

        printf ( "\n" );
        printf ( "QUAD:\n" );
        printf ( "  Estimate the integral of f(x) from A to B.\n" );
        printf ( "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n" );
        printf ( "\n" );
        printf ( "  A        = %f\n", a );
        printf ( "  B        = %f\n", b );
        printf ( "  N        = %d\n", n );
        printf ( "  Exact    = %24.16f\n", exact );
    }


    wtime = omp_get_wtime ( );



//    printf("myrank %d commsize %d\n", myrank, commSize);

    MPI_Bcast( &a, 1, MPI_INT, 0, MPI_COMM_WORLD);

    n_my = (n + commSize - 1)/commSize;
    chunkSize = (b-a)/commSize;
    a_my = a + chunkSize * myrank;
    b_my = a_my + chunkSize;

    total_my = process(a_my,b_my,n_my);
//    printf("my total: %f\n", total_my);
//    printf("my a %f my b %f my n %f\n", a_my, b_my, n_my);

    MPI_Reduce(&total_my, &total, 1, MPI_DOUBLE, MPI_SUM, 0 ,MPI_COMM_WORLD);

    MPI_Finalize();

//    printf("my final total: %f\n", total);

    wtime = omp_get_wtime ( ) - wtime;

    error = fabs ( total - exact );

    if (myrank == 0) {
        printf ( "\n" );
        printf ( "  Estimate = %24.16f\n", total );
        printf ( "  Error    = %e\n", error );
        printf ( "  Time     = %f\n", wtime );
        printf ( "\n" );
        printf ( "  Normal end of execution.\n" );
        printf ( "\n" );
    }

    (*time) = wtime;
    (*result) = total;


    return 0;
}

int sequential ( int argc, char *argv[], double *result, double *time ) {
    double a;
    double b;
    double error;
    double exact = 0.49936338107645674464;
    int i;
    int n;
    double total;
    double wtime;
    double x;

    printf("SEQUENTIAL RUN\n");

    if (argc != 4) {
        n = 10000000;
        a = 0.0;
        b = 10.0;
    } else {
        n = atoi(argv[1]);
        a = atoi(argv[2]);
        b = atoi(argv[3]);
    }

    printf ( "\n" );
    printf ( "QUAD:\n" );
    printf ( "  Estimate the integral of f(x) from A to B.\n" );
    printf ( "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n" );
    printf ( "\n" );
    printf ( "  A        = %f\n", a );
    printf ( "  B        = %f\n", b );
    printf ( "  N        = %d\n", n );
    printf ( "  Exact    = %24.16f\n", exact );

    wtime = omp_get_wtime ( );

    total = 0.0;

    for ( i = 0; i < n; i++ )
    {
        x = ( ( double ) ( n - i - 1 ) * a + ( double ) ( i ) * b ) / ( double ) ( n - 1 );
        total = total + f ( x );
    }

    wtime = omp_get_wtime ( ) - wtime;

    total = ( b - a ) * total / ( double ) n;
    error = fabs ( total - exact );

    printf ( "\n" );
    printf ( "  Estimate = %24.16f\n", total );
    printf ( "  Error    = %e\n", error );
    printf ( "  Time     = %f\n", wtime );
    printf ( "\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );

    (*time) = wtime;
    (*result) = total;

    return 0;
}

double f ( double x ) {
    double pi = 3.141592653589793;
    double value;

    value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

    return value;
}


int main( int argc, char *argv[] ) {

    double sequential_result, parallel_result, sequential_time, parallel_time;


    /* Initialize MPI */

    MPI_Init(&argc, &argv);

    /* Find out my identity in the default communicator */

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (myrank == 0) {
        sequential(argc, argv, &sequential_result, &sequential_time);
    }
    parallel(argc, argv, &parallel_result, &parallel_time);

    if (myrank == 0) {
        compare_and_print(sequential_result, parallel_result, sequential_time, parallel_time);
    }
}

