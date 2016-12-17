# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( int argc, char *argv[] );
double cpu_time ( void );

int main ( int argc, char *argv[] )
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

double cpu_time ( void )
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
