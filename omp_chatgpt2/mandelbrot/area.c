// Correct, but slow
// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <omp.h> 

// #define NPOINTS 2000
// #define MAXITER 2000

// struct complex {
//   double real;
//   double imag;
// };

// int main(int argc, char *argv[]) {
//   int i, j, iter, numoutside = 0, npoints, maxiter;
//   double area, error, ztemp;
//   struct complex z, c;
//   double timer_start, timer_end, elapsed_time;

//   if (argc != 3) {
//     printf("Wrong parameters: npoints and maxiter needed!\n");
//     return 1;
//   }
//   npoints = atoi(argv[1]);
//   maxiter = atoi(argv[2]);  

//   timer_start = omp_get_wtime();

// #pragma omp parallel for default(none) private(i, j, iter, z, c, ztemp) reduction(+:numoutside) shared(npoints, maxiter)
//   for (i = 0; i < npoints; i++) {
//     for (j = 0; j < npoints; j++) {
//       c.real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
//       c.imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
//       z = c;
//       for (iter = 0; iter < maxiter; iter++) {
//         ztemp = (z.real * z.real) - (z.imag * z.imag) + c.real;
//         z.imag = z.real * z.imag * 2 + c.imag; 
//         z.real = ztemp; 
//         if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
//           numoutside++; 
//           break;
//         }
//       }
//     }
//   }

//   timer_end = omp_get_wtime();

//   area = 2.0 * 2.5 * 1.125 * (double)(npoints * npoints - numoutside) / (double)(npoints * npoints);
//   error = area / (double)npoints;

//   printf("Area of Mandelbrot set = %12.8f +/- %12.8f\n", area, error);

//   elapsed_time = timer_end - timer_start;
//   printf("Elapsed time: %.2f\n", elapsed_time);

//   return 0;
// }

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

#define NPOINTS 2000
#define MAXITER 2000

struct complex {
  double real;
  double imag;
};

int main(int argc, char *argv[]) {
  int i, j, iter, numoutside = 0, npoints, maxiter;
  double area, error, ztemp;
  struct complex z, c;
  double timer_start, timer_end, elapsed_time;

  if (argc != 3) {
    printf("Wrong parameters: npoints and maxiter needed!\n");
    return 1;
  }
  npoints = atoi(argv[1]);
  maxiter = atoi(argv[2]);  

  timer_start = omp_get_wtime();

// First answer
//#pragma omp parallel for collapse(2) default(none) private(i, j, iter, z, c, ztemp) reduction(+:numoutside) shared(npoints, maxiter)
// Second answer
#pragma omp parallel for collapse(2) default(none) private(i, j, iter, z, c, ztemp) reduction(+:numoutside) shared(npoints, maxiter) schedule(dynamic)

  for (i = 0; i < npoints; i++) {
    for (j = 0; j < npoints; j++) {
      c.real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
      c.imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
      z = c;
      for (iter = 0; iter < maxiter; iter++) {
        ztemp = (z.real * z.real) - (z.imag * z.imag) + c.real;
        z.imag = z.real * z.imag * 2 + c.imag; 
        z.real = ztemp; 
        if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
          numoutside++; 
          break;
        }
      }
    }
  }

  timer_end = omp_get_wtime();

  area = 2.0 * 2.5 * 1.125 * (double)(npoints * npoints - numoutside) / (double)(npoints * npoints);
  error = area / (double)npoints;

  printf("Area of Mandelbrot set = %12.8f +/- %12.8f\n", area, error);

  elapsed_time = timer_end - timer_start;
  printf("Elapsed time: %.6f\n", elapsed_time);

  return 0;
}
