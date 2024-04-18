// area.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

# define NPOINTS 2000
# define MAXITER 2000


struct complex{
  double real;
  double imag;
};

int main(int argc, char *argv[]){
  int i, j, iter, numoutside = 0, npoints, maxiter;
  double area, error, ztemp;
  struct complex z, c;
  double timer_start, timer_end, elapsed_time;

/*
 *   
 *
 *     Outer loops run over npoints, initialise z=c
 *
 *     Inner loop has the iteration z=z*z+c, and threshold test
 */

  if (argc != 3) {
    printf("Wrong parameters: npoints and maxiter needed!\n");
    return 1;
  }
  npoints = atoi(argv[1]);
  maxiter = atoi(argv[2]);  

  timer_start = omp_get_wtime();

# pragma omp parallel for default(none) shared(npoints, maxiter) \
                           private(i, j, c, z, iter, ztemp) \
                           schedule(dynamic) \
                           reduction(+:numoutside)
  for (i=0; i<npoints; i++) {
    for (j=0; j<npoints; j++) {
      c.real = -2.0+2.5*(double)(i)/(double)(npoints)+1.0e-7;
      c.imag = 1.125*(double)(j)/(double)(npoints)+1.0e-7;
      z=c;
      for (iter=0; iter<maxiter; iter++){
        ztemp=(z.real*z.real)-(z.imag*z.imag)+c.real;
        z.imag=z.real*z.imag*2+c.imag; 
        z.real=ztemp; 
        if ((z.real*z.real+z.imag*z.imag)>4.0e0) {
          numoutside++; 
          break;
        }
      }
    }
  }
  timer_end = omp_get_wtime();

  /*
  *  Calculate area and error and output the results
  */

  area=2.0*2.5*1.125*(double)(npoints*npoints-numoutside)/(double)(npoints*npoints);
  error=area/(double)npoints;

  printf("Area of Mandlebrot set = %12.8f +/- %12.8f\n",area,error);

  elapsed_time = timer_end - timer_start;
  printf("Elapsed time: %.6f\n", elapsed_time);

  return 0;
}
