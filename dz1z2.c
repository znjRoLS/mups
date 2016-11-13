# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

#include "common.h"

int main ( int argc, char *argv[] );
double f ( double x );

int sequential ( int argc, char *argv[], double *result );
int parallel ( int argc, char *argv[], double *result );

double f ( double x ) {
    double pi = 3.141592653589793;
    double value;

    value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

    return value;
}

int sequential ( int argc, char *argv[], double *result ) {
    double a;
    double b;
    double error;
    double exact = 0.49936338107645674464;
    int i;
    int n;
    double total;
    double wtime;
    double x;

    if (argc != 4) {
        n = 100000000;
        a = 0.0;
        b = 10.0;
    } else {
        n = atoi(argv[1]);
        a = atoi(argv[2]);
        b = atoi(argv[3]);
    }

    printf ( "\n" );
    printf ( "QUAD sequential:\n" );
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
    *result = total;
    printf ( "  Error    = %e\n", error );
    printf ( "  Time     = %f\n", wtime );
    printf ( "\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );

    return 0;
}

int parallel ( int argc, char *argv[], double *result ) {
    double a;
    double b;
    double error;
    double exact = 0.49936338107645674464;
    int i;
    int n;
    double total;
    double wtime;
    double x;

    if (argc != 4) {
        n = 100000000;
        a = 0.0;
        b = 10.0;
    } else {
        n = atoi(argv[1]);
        a = atoi(argv[2]);
        b = atoi(argv[3]);
    }

    printf ( "\n" );
    printf ( "QUAD parallel:\n" );
    printf ( "  Estimate the integral of f(x) from A to B.\n" );
    printf ( "  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).\n" );
    printf ( "\n" );
    printf ( "  A        = %f\n", a );
    printf ( "  B        = %f\n", b );
    printf ( "  N        = %d\n", n );
    printf ( "  Exact    = %24.16f\n", exact );

    wtime = omp_get_wtime ( );

    total = 0.0;

#pragma omp parallel for \
    private(x)
    for ( i = 0; i < n; i++ )
    {
        x = ( ( double ) ( n - i - 1 ) * a + ( double ) ( i ) * b ) / ( double ) ( n - 1 );
#pragma omp atomic
        total = total + f ( x );
    }

    wtime = omp_get_wtime ( ) - wtime;

    total = ( b - a ) * total / ( double ) n;
    error = fabs ( total - exact );

    printf ( "\n" );
    printf ( "  Estimate = %24.16f\n", total );
    *result = total;
    printf ( "  Error    = %e\n", error );
    printf ( "  Time     = %f\n", wtime );
    printf ( "\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );

    return 0;
}

int main ( int argc, char *argv[]) {

    double sequential_result, parallel_result;
    int err;

    err = parallel(argc, argv, &parallel_result);
    if (err) { return err; }

    err = sequential(argc, argv, &sequential_result);
    if (err) { return err; }

    finish(sequential_result, parallel_result);
}


