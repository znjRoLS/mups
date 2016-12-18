# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>


#include "common.h"
#include <mpi.h>


int myrank, commSize;

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

    printf("SEQUENTIAL RUN %d\n", myrank);

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
    double my_total_q, my_total_t, my_total_s;
    double wtime_q, wtime_t, wtime_s;
    double x;
    double h;

    double results_values[4], results_times[4];
    double result_value, result_time;

    int n_my;
    double a_my, b_my, chunkSize;

    int myGroupRank, groupCommSize;


    int myColor = (myrank - 1)  % 3;

    printf("PARALLEL RUN %d\n", myrank);

    if (myrank == 0 ) {
        myColor = MPI_UNDEFINED;
    }

    const int roots[4] = {0, 1, 2, 3};

    if (myrank == 0) {
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
    }

// Get the group of processes in MPI_COMM_WORLD
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);


// Construct a group containing all of the prime ranks in world_group
    MPI_Group roots_group;
    MPI_Group_incl(world_group, 4, roots, &roots_group);

// Create a new communicator based on the group
    MPI_Comm roots_comm;
    MPI_Comm_create_group(MPI_COMM_WORLD, roots_group, 0, &roots_comm);


    if (myrank < 4) {
        MPI_Bcast( &a, 1, MPI_DOUBLE, 0, roots_comm);
        MPI_Bcast( &b, 1, MPI_DOUBLE, 0, roots_comm);
        MPI_Bcast( &n, 1, MPI_INT, 0, roots_comm);
    }

//    printf("im %d, and my b is %f\n", myrank, b);

//    if (myrank != 0) {
        MPI_Comm group_comm;

        MPI_Comm_split(
                MPI_COMM_WORLD,
                myColor,
                myrank,
                &group_comm
        );

        if (myColor != MPI_UNDEFINED) {
            MPI_Comm_rank(group_comm, &myGroupRank);
            MPI_Comm_size(group_comm, &groupCommSize);
        }


        switch (myColor) {

            case MPI_UNDEFINED:
                break;

            case 0:
                MPI_Bcast(&a, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&b, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&n, 1, MPI_INT, 0, group_comm);

                // Quadratic rule
                wtime_q = omp_get_wtime ( );

                n_my = (n + groupCommSize - 1)/groupCommSize;
                chunkSize = (b-a)/groupCommSize;
                a_my = a + chunkSize * myGroupRank;
                b_my = a_my + chunkSize;

                my_total_q = quadratic(a_my,b_my,n_my);

//                printf("global rank %d, group rank %d, color %d, my_a %f, my_b %f, my_n %d, my_total_q %f\n", myrank, myGroupRank, myColor, a_my, b_my, n_my, my_total_q);

                MPI_Reduce(&my_total_q, &total_q, 1, MPI_DOUBLE, MPI_SUM, 0 ,group_comm);

                wtime_q = omp_get_wtime ( ) - wtime_q;

                result_value = total_q;
                result_time = wtime_q;

                break;

            case 1:
                MPI_Bcast(&a, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&b, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&n, 1, MPI_INT, 0, group_comm);

                // Trapezoidal rule
                wtime_t = omp_get_wtime ( );

                n_my = (n + groupCommSize - 1)/groupCommSize;
                chunkSize = (b-a)/groupCommSize;
                a_my = a + chunkSize * myGroupRank;
                b_my = a_my + chunkSize;

                my_total_t = trapezoidal(a_my,b_my,n_my);

//                printf("global rank %d, group rank %d, color %d, my_a %f, my_b %f, my_n %d, my_total_t %f\n", myrank, myGroupRank, myColor, a_my, b_my, n_my, my_total_t);

                MPI_Reduce(&my_total_t, &total_t, 1, MPI_DOUBLE, MPI_SUM, 0 ,group_comm);

                wtime_t = omp_get_wtime ( ) - wtime_t;

                result_value = total_t;
                result_time = wtime_t;

                break;

            case 2:
                MPI_Bcast(&a, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&b, 1, MPI_DOUBLE, 0, group_comm);
                MPI_Bcast(&n, 1, MPI_INT, 0, group_comm);

                // Simpson 1/3 rule
                wtime_s = omp_get_wtime ( );

                n_my = (n + groupCommSize - 1)/groupCommSize;
                chunkSize = (b-a)/groupCommSize;
                a_my = a + chunkSize * myGroupRank;
                b_my = a_my + chunkSize;

                my_total_s = simpsons(a_my,b_my,n_my);

//                printf("global rank %d, group rank %d, color %d, my_a %f, my_b %f, my_n %d, my_total_s %f\n", myrank, myGroupRank, myColor, a_my, b_my, n_my, my_total_s);

                MPI_Reduce(&my_total_s, &total_s, 1, MPI_DOUBLE, MPI_SUM, 0 ,group_comm);

                wtime_s = omp_get_wtime ( ) - wtime_s;

                result_value = total_s;
                result_time = wtime_s;

                break;
        }
//    }


//    printf ("yep im herhe %d\n", myrank );
    if (myrank < 4) {

//        printf("yep im herhe2 %d\n", myrank );
//        printf("before gather myrank %d, result_value %f, result_values[0] %f, result_values[1] %f, result_values[2] %f", )

        MPI_Gather(&result_value, 1, MPI_DOUBLE, results_values, 1, MPI_DOUBLE, 0, roots_comm );
        MPI_Gather(&result_time, 1, MPI_DOUBLE, results_times, 1, MPI_DOUBLE, 0, roots_comm );
    }
//    printf("yep im herhe3 %d\n", myrank );

    for (int i = 0 ; i < 3; i ++) {
        results[i].value = results_values[i+1];
        results[i].time = results_times[i+1];
    }

    if (myrank == 0) {
        printf ( "\n" );
        printf ( "  Estimate quadratic rule = %24.16f\n", results[0].value );
        printf ( "  Estimate trapezoidal rule = %24.16f\n", results[1].value );
        printf ( "  Estimate Simpson 1/3 rule = %24.16f\n", results[2].value );
        printf ( "  Time quadratic rule = %f\n", results[0].time );
        printf ( "  Time trapezoidal rule = %f\n", results[1].time );
        printf ( "  Time Simpson 1/3 rule = %f\n", results[2].time );
        printf ( "\n" );
        printf ( "  Normal end of execution.\n" );
        printf ( "\n" );
    }

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

    /* Initialize MPI */

    MPI_Init(&argc, &argv);

    /* Find out my identity in the default communicator */

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (myrank == 0) {
        sequential(argc, argv, seq_results);
    }

    parallel(argc, argv, par_results);

//    printf("izaso sam %d\n", myrank);

    if (myrank == 0) {
        compare_and_print(seq_results[0], par_results[0], "quadratic");
        compare_and_print(seq_results[1], par_results[1], "trapezoidal");
        compare_and_print(seq_results[2], par_results[2], "simpsons");
    }

    MPI_Finalize();
}

