# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

#include "common.h"
#include <string.h>
#include <mpi.h>

int main ( int argc, char *argv[] );
double cpu_time ( void );

int myrank, commSize;

double cpu_time ( void )
{
    double value;

    value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

    return value;
}

int sequential ( int argc, char *argv[], Result_Vect *result )
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

    if (argc != 5) {
        printf("Wrong number of arguments!\n");
        return 1;
    } else {
        success = sscanf ( argv[1], "%d", &M );
        success += sscanf ( argv[2], "%d", &N );
        success += sscanf ( argv[3], "%lf", &epsilon );
        success += sscanf ( argv[4], "%s", output_file );

        if (success != 4) {
            printf("Wrong arguments!\n");
            return 2;
        }
    }

    printf("SEQUENTIAL RUN %d\n", myrank);
    printf ( "\n" );
    printf ( "HEATED_PLATE\n" );
    printf ( "  C version\n" );
    printf ( "  A program to solve for the steady state temperature distribution\n" );
    printf ( "  over a rectangular plate.\n" );
    printf ( "\n" );
    printf ( "  Spatial grid of %d by %d points.\n", M, N );
    printf ( "\n" );
    printf ( "  The iteration will be repeated until the change is <= %f\n", epsilon );
    diff = epsilon;
    printf ( "\n" );
    printf ( "  The steady state solution will be written to %s.\n", output_file );

    u = (double **) malloc(M * sizeof(double*));
    for (i = 0; i < M; i++)
        u[i] = (double *) malloc(N * sizeof(double));

    w = (double **) malloc(M * sizeof(double*));
    for (i = 0; i < M; i++)
        w[i] = (double *) malloc(N * sizeof(double));

/*
  Set the boundary values, which don't change.
*/
    for ( i = 1; i < M - 1; i++ )
    {
        w[i][0] = 100.0;
    }
    for ( i = 1; i < M - 1; i++ )
    {
        w[i][N-1] = 100.0;
    }
    for ( j = 0; j < N; j++ )
    {
        w[M-1][j] = 100.0;
    }
    for ( j = 0; j < N; j++ )
    {
        w[0][j] = 0.0;
    }
/*
  Average the boundary values, to come up with a reasonable
  initial value for the interior.
*/
    mean = 0.0;
    for ( i = 1; i < M - 1; i++ )
    {
        mean = mean + w[i][0];
    }
    for ( i = 1; i < M - 1; i++ )
    {
        mean = mean + w[i][N-1];
    }
    for ( j = 0; j < N; j++ )
    {
        mean = mean + w[M-1][j];
    }
    for ( j = 0; j < N; j++ )
    {
        mean = mean + w[0][j];
    }
    mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
/*
  Initialize the interior solution to the mean value.
*/
    for ( i = 1; i < M - 1; i++ )
    {
        for ( j = 1; j < N - 1; j++ )
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
    printf ( "\n" );
    printf ( " Iteration  Change\n" );
    printf ( "\n" );
    ctime1 = cpu_time ( );

    while ( epsilon <= diff )
    {
/*
  Save the old solution in U.
*/
        for ( i = 0; i < M; i++ )
        {
            for ( j = 0; j < N; j++ )
            {
                u[i][j] = w[i][j];
            }
        }
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
        diff = 0.0;
        for ( i = 1; i < M - 1; i++ )
        {
            for ( j = 1; j < N - 1; j++ )
            {
                w[i][j] = ( u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1] ) / 4.0;

                if ( diff < fabs ( w[i][j] - u[i][j] ) )
                {
                    diff = fabs ( w[i][j] - u[i][j] );
                }
            }
        }
        iterations++;
        if ( iterations == iterations_print )
        {
            printf ( "  %8d  %f\n", iterations, diff );
            iterations_print = 2 * iterations_print;
        }
    }
    ctime2 = cpu_time ( );
    ctime = ctime2 - ctime1;

    printf ( "\n" );
    printf ( "  %8d  %f\n", iterations, diff );
    printf ( "\n" );
    printf ( "  Error tolerance achieved.\n" );
    printf ( "  CPU time = %f\n", ctime );

    result->time = ctime;

    result->val_size = M * N;
    result->value = (double*) malloc (sizeof(double) * result->val_size);

/*
  Write the solution to the output file.
*/
    fp = fopen ( output_file, "w" );

    fprintf ( fp, "%d\n", M );
    fprintf ( fp, "%d\n", N );

    for ( i = 0; i < M; i++ )
    {
        for ( j = 0; j < N; j++)
        {
            fprintf ( fp, "%6.2f ", w[i][j] );
            result->value[i * N + j] = w[i][j];
        }
        fputc ( '\n', fp);
    }
    fclose ( fp );

    printf ( "\n" );
    printf ("  Solution written to the output file %s\n", output_file );
/*
  All done!
*/
    printf ( "\n" );
    printf ( "HEATED_PLATE:\n" );
    printf ( "  Normal end of execution.\n" );

    return 0;

}

int parallel ( int argc, char *argv[], Result_Vect *result )
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
    char output_file_par[80];
    int success;

    double **u;
    double **w;

    if (myrank == 0) {

        if (argc != 5) {
            printf("Wrong number of arguments!\n");
            return 1;
        } else {
            success = sscanf ( argv[1], "%d", &M );
            success += sscanf ( argv[2], "%d", &N );
            success += sscanf ( argv[3], "%lf", &epsilon );
            success += sscanf ( argv[4], "%s", output_file );

            for (int i = strlen(output_file)-1; i >= 0; i--) {
                output_file_par[i] = output_file[i];
                output_file_par[0] = 'p';
            }

            if (success != 4) {
                printf("Wrong arguments!\n");
                return 2;
            }
        }

        printf("PARALLEL RUN %d\n", myrank);

        printf ( "\n" );
        printf ( "HEATED_PLATE\n" );
        printf ( "  C version\n" );
        printf ( "  A program to solve for the steady state temperature distribution\n" );
        printf ( "  over a rectangular plate.\n" );
        printf ( "\n" );
        printf ( "  Spatial grid of %d by %d points.\n", M, N );
        printf ( "\n" );
        printf ( "  The iteration will be repeated until the change is <= %f\n", epsilon );
        diff = epsilon;
        printf ( "\n" );
        printf ( "  The steady state solution will be written to %s.\n", output_file );

        u = (double **) malloc(M * sizeof(double*));
        for (i = 0; i < M; i++)
            u[i] = (double *) malloc(N * sizeof(double));

        w = (double **) malloc(M * sizeof(double*));
        for (i = 0; i < M; i++)
            w[i] = (double *) malloc(N * sizeof(double));

    /*
      Set the boundary values, which don't change.
    */
        for ( i = 1; i < M - 1; i++ )
        {
            w[i][0] = 100.0;
        }
        for ( i = 1; i < M - 1; i++ )
        {
            w[i][N-1] = 100.0;
        }
        for ( j = 0; j < N; j++ )
        {
            w[M-1][j] = 100.0;
        }
        for ( j = 0; j < N; j++ )
        {
            w[0][j] = 0.0;
        }
    /*
      Average the boundary values, to come up with a reasonable
      initial value for the interior.
    */
        mean = 0.0;
        for ( i = 1; i < M - 1; i++ )
        {
            mean = mean + w[i][0];
        }
        for ( i = 1; i < M - 1; i++ )
        {
            mean = mean + w[i][N-1];
        }
        for ( j = 0; j < N; j++ )
        {
            mean = mean + w[M-1][j];
        }
        for ( j = 0; j < N; j++ )
        {
            mean = mean + w[0][j];
        }
        mean = mean / ( double ) ( 2 * M + 2 * N - 4 );
    /*
      Initialize the interior solution to the mean value.
    */
        for ( i = 1; i < M - 1; i++ )
        {
            for ( j = 1; j < N - 1; j++ )
            {
                w[i][j] = mean;
            }
        }

    }
/*
  iterate until the  new solution W differs from the old solution U
  by no more than EPSILON.
*/

    MPI_Datatype rowType;
    //MPI_Type_contiguous(N, MPI_DOUBLE, &rowType);
    MPI_Type_contiguous(500, MPI_DOUBLE, &rowType);
//    printf("create type with size %d\n", N);
    MPI_Type_commit(&rowType);

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int columnsNum = M-2;

    int chunkSize = (columnsNum + commSize - 1) / commSize;

    int myStart = myrank * chunkSize + 1;
    int myEnd = myStart + chunkSize;
    if (myEnd > M-1)
        myEnd = M-1;

    int myStartWide = myStart - 1;
    int myEndWide = myEnd + 1;

    int myNumColumns = myEndWide - myStartWide;

    double** myColumns = (double**) malloc (sizeof(double*) * myNumColumns);
    for (int i = 0; i < myEndWide - myStartWide; i ++) {
        myColumns[i] = (double*) malloc(sizeof(double) * N);
    }
//    printf("Creted columns with size %d\n", N);

    //transfer stuff
    if (myrank == 0){
        for (int i = 0; i < myEndWide - myStartWide; i ++) {
            for (int  j = 0 ; j < N ; j ++) {
                myColumns[i][j] = w[i][j];
            }
        }


        for (int irank = 1; irank < commSize; irank ++) {
            int otherStart = irank * chunkSize + 1;
            int otherEnd = otherStart + chunkSize;
            if (otherEnd > M-1)
                otherEnd = M-1;

            int otherStartWide = otherStart - 1;
            int otherEndWide = otherEnd + 1;


//            printf("irank %d, otherStartW %d, otherEndWide %d\n", irank, otherStartWide, otherEndWide);

            for (int i = otherStartWide; i < otherEndWide; i ++) {
//                printf("Sending column %d to %d\n", i, irank);
                MPI_Send(w[i], 1, rowType, irank, 123, MPI_COMM_WORLD);
            }
        }
    }
    else {
        for (int i = 0; i < myEndWide - myStartWide; i ++) {
//            printf("Receiveing column %d, im %d\n", i, myrank);
            MPI_Status status;
            MPI_Recv(myColumns[i], 1, rowType, 0, 123, MPI_COMM_WORLD, &status);
        }
    }

    iterations = 0;
    iterations_print = 1;
    printf ( "\n" );
    printf ( " Iteration  Change\n" );
    printf ( "\n" );
    ctime1 = cpu_time ( );

    double** myOldColumns = (double**) malloc (sizeof(double*) * myNumColumns);
    for (int i = 0; i < myNumColumns; i ++) {
        myOldColumns[i] = (double*) malloc(sizeof(double) * N);
    }

    double total_diff = epsilon;

    while ( epsilon <= total_diff )
    {
//        printf("jel sam opet uso?\n");
/*
  Save the old solution in U.
*/
        for ( i = 0; i < myNumColumns; i++ )
        {
            for ( j = 0; j < N; j++ )
            {
                myOldColumns[i][j] = myColumns[i][j];
            }
        }

        if (myrank == 0) {
            //communicate only with right one
            MPI_Status status;

            MPI_Send(myOldColumns[myNumColumns-2], 1, rowType, 1, 0, MPI_COMM_WORLD);
            MPI_Recv(myOldColumns[myNumColumns-1], 1, rowType, 1, 0, MPI_COMM_WORLD, &status);

        }
        else if (myrank == commSize - 1) {
            //communicate only with left one
            MPI_Status status;

            MPI_Send(myOldColumns[1], 1, rowType, myrank-1, 0, MPI_COMM_WORLD);
            MPI_Recv(myOldColumns[0], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &status);
        }

        else {
            //communicate with both
            MPI_Status status;

            MPI_Send(myOldColumns[1], 1, rowType, myrank-1, 0, MPI_COMM_WORLD);
            MPI_Send(myOldColumns[myNumColumns-2], 1, rowType, myrank+1, 0, MPI_COMM_WORLD);
            MPI_Recv(myOldColumns[0], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(myOldColumns[myNumColumns-1], 1, rowType, myrank+1, 0, MPI_COMM_WORLD, &status);

        }
/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
        diff = 0.0;
        total_diff = 0.0;
        for ( i = 1; i < myNumColumns - 1; i++ )
        {
            for ( j = 1; j < N - 1; j++ )
            {
                myColumns[i][j] = ( myOldColumns[i-1][j] + myOldColumns[i+1][j] + myOldColumns[i][j-1] + myOldColumns[i][j+1] ) / 4.0;

                if ( diff < fabs ( myColumns[i][j] - myOldColumns[i][j] ) )
                {
                    diff = fabs ( myColumns[i][j] - myOldColumns[i][j] );
                }
            }
        }
        iterations++;

        MPI_Allreduce(&diff, &total_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//        total_diff /= commSize;
//
//        MPI_Bcast(&total_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//        printf("herhe %d - %d - %f\n", myrank, iterations);

        if ( iterations == iterations_print )
        {
//            printf ( "  %8d  %f\n", iterations, total_diff );
            iterations_print = 2 * iterations_print;
        }
    }
    ctime2 = cpu_time ( );
    ctime = ctime2 - ctime1;

//    printf("kk?\n");

    //transfer stuff back
    if (myrank == 0){
//        printf("tuj sam\n");
        for (int i = 0; i < myEndWide - myStartWide; i ++) {
            for (int  j = 0 ; j < N ; j ++) {
                w[i][j] = myColumns[i][j];
            }
        }


        for (int irank = 1; irank < commSize; irank ++) {
            int otherStart = irank * chunkSize + 1;
            int otherEnd = otherStart + chunkSize;
            if (otherEnd > M-1)
                otherEnd = M-1;

            int otherStartWide = otherStart - 1;
            int otherEndWide = otherEnd + 1;


//            printf("irank %d, otherStartW %d, otherEndWide %d\n", irank, otherStartWide, otherEndWide);

            for (int i = otherStartWide+1; i < otherEndWide; i ++) {
//                printf("Sending column %d to %d\n", i, irank);
                MPI_Status status;
                MPI_Recv(w[i], 1, rowType, irank, 123, MPI_COMM_WORLD, &status);
            }
        }
    }
    else {
//        printf("a i ja tuj sam\n");
        for (int i = 1; i < myEndWide - myStartWide; i ++) {
//            printf("Receiveing column %d, im %d\n", i, myrank);
            MPI_Send(myColumns[i], 1, rowType, 0, 123, MPI_COMM_WORLD);
        }
    }

    if (myrank == 0) {
        printf ( "\n" );
        printf ( "  %8d  %f\n", iterations, diff );
        printf ( "\n" );
        printf ( "  Error tolerance achieved.\n" );
        printf ( "  CPU time = %f\n", ctime );

        result->time = ctime;

        result->val_size = M * N;
        result->value = (double*) malloc (sizeof(double) * result->val_size);
/*
  Write the solution to the output file.
*/
        fp = fopen ( output_file_par, "w" );

        fprintf ( fp, "%d\n", M );
        fprintf ( fp, "%d\n", N );

        for ( i = 0; i < M; i++ )
        {
            for ( j = 0; j < N; j++)
            {
                fprintf ( fp, "%6.2f ", w[i][j] );
                result->value[i * N + j] = w[i][j];
            }
            fputc ( '\n', fp);
        }
        fclose ( fp );

        printf ( "\n" );
        printf ("  Solution written to the output file %s\n", output_file_par );
/*
  All done!
*/
        printf ( "\n" );
        printf ( "HEATED_PLATE:\n" );
        printf ( "  Normal end of execution.\n" );
    }

    return 0;
}

int main (int argc, char *argv[]) {
    Result_Vect seq_result;
    Result_Vect par_result;

    /* Initialize MPI */

    MPI_Init(&argc, &argv);

    /* Find out my identity in the default communicator */

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    if (myrank == 0) {
        sequential(argc, argv, &seq_result);
    }

    parallel(argc, argv, &par_result);

    if (myrank == 0) {
        compare_and_print_vect(seq_result, par_result, "HEATED PLATE TEST");
    }

    MPI_Finalize();
}

//int main (int argc, char *argv[]) {
//    Result_Vect seq_result;
//    Result_Vect par_result;
//
//    /* Initialize MPI */
//
//    MPI_Init(&argc, &argv);
//
//    /* Find out my identity in the default communicator */
//
//    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
//
//    double **a;
//    a = (double **) malloc(sizeof(double *) * 2);
//    a[0] = (double *) malloc(sizeof(double) * 500);
//
//    double **b;
//    b = (double **) malloc(sizeof(double *) * 2);
//    b[0] = (double *) malloc(sizeof(double) * 500);
//
//    if (myrank == 0) {
//        for (int i = 0; i < 500; i++)
//            a[0][i] = i;
//    }
//
//    MPI_Datatype rowType;
//
//    MPI_Type_contiguous(500, MPI_DOUBLE, &rowType);
//    MPI_Type_commit(&rowType);
//
//    if (myrank == 0) {
//        MPI_Send(a[0], 1, rowType, 1, 123, MPI_COMM_WORLD);
//        printf("snettt");
//    }
//    if (myrank == 1 ) {
//        MPI_Status status;
//        MPI_Recv(b[0], 1, rowType, 0, 123, MPI_COMM_WORLD, &status);
//        printf("receiveddd");
//    }
//
//    MPI_Finalize();
//}

