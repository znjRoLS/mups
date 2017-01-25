# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include <omp.h>

int main ( int argc, char *argv[] );
double f ( double x );

int main ( int argc, char *argv[] ) {
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

  total_q = 0.0;

  for ( i = 0; i < n; i++ )
  {
    x = ( ( double ) ( n - i - 1 ) * a + ( double ) ( i ) * b ) / ( double ) ( n - 1 );
    total_q = total_q + f ( x );
  }

  wtime_q = omp_get_wtime ( ) - wtime_q;

  total_q = ( b - a ) * total_q / ( double ) n;

  
  // Trapezoidal rule  
  h = (b - a) / n;
  
  wtime_t = omp_get_wtime ( );

  total_t = 0.0;
      
  for ( i = 0; i < n; i++ )
  {
    x = a + i * h;
    if (i > 0 && i < n - 1) 
      total_t = total_t + f( x );
    else 
      total_t = total_t + 0.5 * f( x );
  }

  total_t = h * total_t;

  wtime_t = omp_get_wtime ( ) - wtime_t;

  // Simpson 1/3 rule  

  h = (b - a) / n;
  
  wtime_s = omp_get_wtime ( );

  total_s = 0.0;
      
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

  wtime_s = omp_get_wtime ( ) - wtime_s;

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

