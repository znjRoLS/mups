# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>


#include "common.h"
#include <mpi.h>

int main ( int argc, char *argv[] );
double f ( double x );

double quadratic (double a, double b, int n){
    int i;
    double x;
    double total_q = 0.0;

    for ( i = 0; i < n; i++ )
    {
        x = ( ( double ) ( n - i - 1 ) * a + ( double ) ( i ) * b ) / ( double ) ( n - 1 );
        total_q = total_q + f ( x );
    }

    total_q = ( b - a ) * total_q / ( double ) n;

    return total_q;
}

double trapezoidal (double a, double b, int n) {
    double total_t = 0.0;
    int i;
    double x;
    double h = (b - a) / n;

    for ( i = 0; i < n; i++ )
    {
        x = a + i * h;
        if (i > 0 && i < n - 1)
            total_t = total_t + f( x );
        else
            total_t = total_t + 0.5 * f( x );
    }

    total_t = h * total_t;

    return total_t;
}

double simpsons(double a, double b, int n) {
    double h = (b - a) / n;
    double total_s = 0.0;
    int i;
    double x;

    for ( i = 0; i < n; i++ )
    {
        x = a + i * h;
        if (i == 0 || i == n - 1)
            total_s = total_s + f( x );
        else if (i % 2 == 1)
            total_s = total_s + 4 * f( x );
        else
            total_s = total_s + 2 * f( x );
    }

    total_s = h / 3 * total_s;

    return total_s;
}

int sequential ( int argc, char *argv[], Result *results ) {
    double a;
    double b;
    double error;
    int i;
    int n;
    double total_q, total_t, total_s;
    double wtime_q, wtime_t, wtime_s;
    double x;
    double h;

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


    // Quadratic rule
    wtime_q = omp_get_wtime ( );

    total_q = quadratic(a,b,n);

    wtime_q = omp_get_wtime ( ) - wtime_q;

    results[0].value = total_q;
    results[0].time = wtime_q;

    // Trapezoidal rule
    wtime_t = omp_get_wtime ( );

    total_t = trapezoidal(a, b, n);

    wtime_t = omp_get_wtime ( ) - wtime_t;

    results[1].value = total_t;
    results[1].time = wtime_t;

    // Simpson 1/3 rule
    wtime_s = omp_get_wtime ( );

    total_s = simpsons(a,b,n);

    wtime_s = omp_get_wtime ( ) - wtime_s;

    results[2].value = total_s;
    results[2].time = wtime_s;

    printf ( "\n" );
    printf ( "  Estimate quadratic rule = %24.16f\n", total_q );
    printf ( "  Estimate trapezoidal rule = %24.16f\n", total_t );
    printf ( "  Estimate Simpson 1/3 rule = %24.16f\n", total_s );
    printf ( "  Time quadratic rule = %f\n", wtime_q );
    printf ( "  Time trapezoidal rule = %f\n", wtime_t );
    printf ( "  Time Simpson 1/3 rule = %f\n", wtime_s );
    printf ( "\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );

    return 0;
}

int parallel ( int argc, char *argv[],  Result *results ) {
    double a;
    double b;
    double error;
    int i;
    int n;
    double total_q, total_t, total_s;
    double wtime_q, wtime_t, wtime_s;
    double x;
    double h;

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


    // Quadratic rule
    wtime_q = omp_get_wtime ( );

    total_q = quadratic(a,b,n);

    wtime_q = omp_get_wtime ( ) - wtime_q;

    results[0].value = total_q;
    results[0].time = wtime_q;

    // Trapezoidal rule
    wtime_t = omp_get_wtime ( );

    total_t = trapezoidal(a, b, n);

    wtime_t = omp_get_wtime ( ) - wtime_t;

    results[1].value = total_t;
    results[1].time = wtime_t;

    // Simpson 1/3 rule
    wtime_s = omp_get_wtime ( );

    total_s = simpsons(a,b,n);

    wtime_s = omp_get_wtime ( ) - wtime_s;

    results[2].value = total_s;
    results[2].time = wtime_s;

    printf ( "\n" );
    printf ( "  Estimate quadratic rule = %24.16f\n", total_q );
    printf ( "  Estimate trapezoidal rule = %24.16f\n", total_t );
    printf ( "  Estimate Simpson 1/3 rule = %24.16f\n", total_s );
    printf ( "  Time quadratic rule = %f\n", wtime_q );
    printf ( "  Time trapezoidal rule = %f\n", wtime_t );
    printf ( "  Time Simpson 1/3 rule = %f\n", wtime_s );
    printf ( "\n" );
    printf ( "  Normal end of execution.\n" );
    printf ( "\n" );

    return 0;
}

double f ( double x ) {
    double pi = 3.141592653589793;
    double value;

    value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

    return value;
}

int main (int argc, char *argv[]) {

    Result seq_results[3];
    Result par_results[3];

    sequential(argc, argv, seq_results);
    parallel(argc, argv, par_results);

    compare_and_print_result(seq_results[0], par_results[0], "quadratic");
    compare_and_print_result(seq_results[1], par_results[1], "trapezoidal");
    compare_and_print_result(seq_results[2], par_results[2], "simpsons");

}

