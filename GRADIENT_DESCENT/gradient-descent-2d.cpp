/*
Gradient descent example for a 2-D function. Minima of f(x,y) = 4x(x-y) + 2y^2. Located at (x,y)=(0,0)
compile as
g++ gradient-descent-2d.cpp -o gradient-descent-2d
run as
gradient-descent-2d <starting x value> <starting y value>
*/

#include <stdio.h>
#include <stdlib.h>

double funct1(double x, double y);

struct gradient
{
double xval;
double yval;
} ;

gradient funct1gradient(double x, double y);

int main(int argc, char *argv[])
{
  double x = atof(argv[1]);
  double y = atof(argv[2]);

  printf("intput args [ x = %f], [y = %f]\n", x, y);

  double h = 0.1;

  for (int idx=1; idx <= 100; idx++)
  {
      gradient g = funct1gradient(x, y);
      x = x - h * g.xval;
      y = y - h * g.yval;
  }

  printf("final point [ x = %f], [y = %f]\n", x, y);
  printf("function val = %f\n", funct1(x, y));
}

double funct1(double x, double y)
{
    return 4*x*(x-y) + 2*y*y ;
}

gradient funct1gradient(double x, double y)
{
    gradient retval = { 0.0, 0.0 };
    retval.xval = 4*( 2*x - y );
    retval.yval = 4*( y - x );

    return retval;
}


