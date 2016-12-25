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

            for (i = strlen(output_file)-1; i >= 0; i--) {
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
    int columnsNum, chunkSize, myStart, myEnd, myStartWide, myEndWide, myNumColumns;
    double ** myColumns, ** myOldColumns;
    double total_diff;
    MPI_Datatype rowType;


    //MPI_Type_contiguous(N, MPI_DOUBLE, &rowType);
    MPI_Type_contiguous(500, MPI_DOUBLE, &rowType);
//    printf("create type with size %d\n", N);
    MPI_Type_commit(&rowType);

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    columnsNum = M-2;

    chunkSize = (columnsNum + commSize - 1) / commSize;

    myStart = myrank * chunkSize + 1;
    myEnd = myStart + chunkSize;
    if (myEnd > M-1)
        myEnd = M-1;

    myStartWide = myStart - 1;
    myEndWide = myEnd + 1;

    myNumColumns = myEndWide - myStartWide;

    myColumns = (double**) malloc (sizeof(double*) * myNumColumns);
    for (i = 0; i < myEndWide - myStartWide; i ++) {
        myColumns[i] = (double*) malloc(sizeof(double) * N);
    }
//    printf("Creted columns with size %d\n", N);

    //transfer stuff
    if (myrank == 0){
        int irank;
        for (i = 0; i < myEndWide - myStartWide; i ++) {
            for ( j = 0 ; j < N ; j ++) {
                myColumns[i][j] = w[i][j];
            }
        }


        for ( irank = 1; irank < commSize; irank ++) {
            int otherStart, otherEnd, otherStartWide, otherEndWide, nekoi;
            otherStart = irank * chunkSize + 1;
            otherEnd = otherStart + chunkSize;
            if (otherEnd > M-1)
                otherEnd = M-1;

            otherStartWide = otherStart - 1;
            otherEndWide = otherEnd + 1;


//            printf("irank %d, otherStartW %d, otherEndWide %d\n", irank, otherStartWide, otherEndWide);

            for (nekoi = otherStartWide; nekoi < otherEndWide; nekoi ++) {
//                printf("Sending column %d to %d\n", i, irank);
                MPI_Send(w[nekoi], 1, rowType, irank, 123, MPI_COMM_WORLD);
            }
        }
    }
    else {
        for (i = 0; i < myEndWide - myStartWide; i ++) {
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

    myOldColumns = (double**) malloc (sizeof(double*) * myNumColumns);
    for (i = 0; i < myNumColumns; i ++) {
        myOldColumns[i] = (double*) malloc(sizeof(double) * N);
    }

    total_diff = epsilon;

    while ( epsilon <= total_diff )
    {
        MPI_Request requestsend1, requestsend2, request1, request2;
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

            MPI_Isend(myOldColumns[myNumColumns-2], 1, rowType, 1, 0, MPI_COMM_WORLD, &requestsend1);
            MPI_Irecv(myOldColumns[myNumColumns-1], 1, rowType, 1, 0, MPI_COMM_WORLD, &request1);

        }
        else if (myrank == commSize - 1) {
            //communicate only with left one

            MPI_Isend(myOldColumns[1], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &requestsend1);
            MPI_Irecv(myOldColumns[0], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &request1);
        }

        else {
            //communicate with both

            MPI_Isend(myOldColumns[1], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &requestsend1);
            MPI_Isend(myOldColumns[myNumColumns-2], 1, rowType, myrank+1, 0, MPI_COMM_WORLD, &requestsend2);
            MPI_Irecv(myOldColumns[0], 1, rowType, myrank-1, 0, MPI_COMM_WORLD, &request1);
            MPI_Irecv(myOldColumns[myNumColumns-1], 1, rowType, myrank+1, 0, MPI_COMM_WORLD, &request2);

        }

/*
  Determine the new estimate of the solution at the interior points.
  The new solution W is the average of north, south, east and west neighbors.
*/
        diff = 0.0;
        total_diff = 0.0;
        for ( i = 2; i < myNumColumns - 2; i++ )
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

        if (myrank == 0 || myrank == commSize - 1) {
            MPI_Status status;

            MPI_Wait(&request1, &status);
        }
        else {
            MPI_Status status;

            MPI_Wait(&request1, &status);
            MPI_Wait(&request2, &status);
        }

        for ( j = 1; j < N - 1; j++ )
        {
            myColumns[1][j] = ( myOldColumns[0][j] + myOldColumns[2][j] + myOldColumns[1][j-1] + myOldColumns[1][j+1] ) / 4.0;

            if ( diff < fabs ( myColumns[1][j] - myOldColumns[1][j] ) )
            {
                diff = fabs ( myColumns[1][j] - myOldColumns[1][j] );
            }
        }

        for ( j = 1; j < N - 1; j++ )
        {
            myColumns[myNumColumns-2][j] = ( myOldColumns[myNumColumns-3][j] + myOldColumns[myNumColumns-1][j] + myOldColumns[myNumColumns-2][j-1] + myOldColumns[myNumColumns-2][j+1] ) / 4.0;

            if ( diff < fabs ( myColumns[myNumColumns-2][j] - myOldColumns[myNumColumns-2][j] ) )
            {
                diff = fabs ( myColumns[myNumColumns-2][j] - myOldColumns[myNumColumns-2][j] );
            }
        }


        iterations++;

        MPI_Allreduce(&diff, &total_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//        total_diff /= commSize;
//
//        MPI_Bcast(&total_diff, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

//        printf("herhe %d - %d - %f\n", myrank, iterations);

        if (myrank == 0 || myrank == commSize - 1) {
            MPI_Status status;

            MPI_Wait(&requestsend1, &status);
        }
        else {
            MPI_Status status;

            MPI_Wait(&requestsend1, &status);
            MPI_Wait(&requestsend2, &status);

        }

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
        int irank;
//        printf("tuj sam\n");
        for (i = 0; i < myEndWide - myStartWide; i ++) {
            for ( j = 0 ; j < N ; j ++) {
                w[i][j] = myColumns[i][j];
            }
        }


        for (irank = 1; irank < commSize; irank ++) {
            int otherStart, otherEnd, otherStartWide, otherEndWide, nekoi;
            otherStart = irank * chunkSize + 1;
            otherEnd = otherStart + chunkSize;
            if (otherEnd > M-1)
                otherEnd = M-1;

            otherStartWide = otherStart - 1;
            otherEndWide = otherEnd + 1;


//            printf("irank %d, otherStartW %d, otherEndWide %d\n", irank, otherStartWide, otherEndWide);

            for (nekoi = otherStartWide+1; nekoi < otherEndWide; nekoi ++) {
//                printf("Sending column %d to %d\n", i, irank);
                MPI_Status status;
                MPI_Recv(w[nekoi], 1, rowType, irank, 123, MPI_COMM_WORLD, &status);
            }
        }
    }
    else {
//        printf("a i ja tuj sam\n");
        for (i = 1; i < myEndWide - myStartWide; i ++) {
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

