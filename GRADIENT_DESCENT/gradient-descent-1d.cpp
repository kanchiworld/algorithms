/*
Gradient descent example. Minima of a parabola f(x) = x^2 - 4x + 1. Located at x=2.0
compile as
g++ gradient-descent-1d.cpp -o gradient-descent-1d
run as
gradient-descent-1d <starting x value>
*/

#include <stdio.h>
#include <stdlib.h>
double funct1(double x );
double funct1gradient(double x );

int main(int argc, char *argv[])
{
  double x = atof(argv[1]);
  double h = 0.1;
  printf("starting point [ x = %f] \n", x);
  
  for (int idx=1; idx <= 100; idx++)
  {
      x = x - h * funct1gradient(x);
      printf("step number %d => next x = %f\n", idx, x);
  }
}

double funct1(double x )
{
    return x*x - 4*x + 1 ;
}

double funct1gradient(double x )
{
    return 2*x - 4;
}


