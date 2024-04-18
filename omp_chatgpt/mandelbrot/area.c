// area.c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h> 

//# define NPOINTS 2000
//# define MAXITER 2000


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
  // for (i=0; i<npoints; i++) {
  //   for (j=0; j<npoints; j++) {
  //     c.real = -2.0+2.5*(double)(i)/(double)(npoints)+1.0e-7;
  //     c.imag = 1.125*(double)(j)/(double)(npoints)+1.0e-7;
  //     z=c;
  //     for (iter=0; iter<maxiter; iter++){
  //       ztemp=(z.real*z.real)-(z.imag*z.imag)+c.real;
  //       z.imag=z.real*z.imag*2+c.imag; 
  //       z.real=ztemp; 
  //       if ((z.real*z.real+z.imag*z.imag)>4.0e0) {
  //         numoutside++; 
  //         break;
  //       }
  //     }
  //   }
  // }

  // Does not compile due to ztem mistaken for complex number
  //   // Parallelize the outer loop
  // #pragma omp parallel for private(z, c, ztemp) reduction(+:numoutside)
  // for (int i = 0; i < npoints; i++) {
  //     for (int j = 0; j < npoints; j++) {
  //         // Compute complex numbers c and z for each point
  //         c.real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
  //         c.imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
  //         z = c;

  //         // Iterate to determine if the point is in the Mandelbrot set
  //         for (int iter = 0; iter < maxiter; iter++) {
  //             ztemp.real = (z.real * z.real) - (z.imag * z.imag) + c.real;
  //             ztemp.imag = z.real * z.imag * 2 + c.imag;
  //             z = ztemp;

  //             // Check if the point is outside the Mandelbrot set
  //             if ((z.real * z.real + z.imag * z.imag) > 4.0e0) {
  //                 #pragma omp atomic
  //                 numoutside++;
  //                 break;
  //             }
  //         }
  //     }
  // }

    // double ztemp_real, ztemp_imag;
    // double c_real, c_imag;

    // // Parallelize the outer two loops
    // #pragma omp parallel for private(c_real, c_imag, ztemp_real, ztemp_imag) reduction(+:numoutside)
    // for (int i = 0; i < npoints; i++) {
    //     for (int j = 0; j < npoints; j++) {
    //         // Define c for each point
    //         c_real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
    //         c_imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
            
    //         // Initialize z for each point
    //         double z_real = c_real;
    //         double z_imag = c_imag;

    //         // Perform iterations for each point
    //         for (int iter = 0; iter < maxiter; iter++) {
    //             ztemp_real = (z_real * z_real) - (z_imag * z_imag) + c_real;
    //             ztemp_imag = z_real * z_imag * 2 + c_imag;
    //             z_real = ztemp_real;
    //             z_imag = ztemp_imag;
                
    //             // Check if z is outside the boundary
    //             if ((z_real * z_real + z_imag * z_imag) > 4.0e0) {
    //                 #pragma omp atomic
    //                 numoutside++;
    //                 break;
    //             }
    //         }
    //     }
    // }

    // double ztemp_real, ztemp_imag;
    // double c_real, c_imag;

    // // Parallelize the outer two loops and use reduction for numoutside
    // #pragma omp parallel for private(c_real, c_imag, ztemp_real, ztemp_imag) reduction(+:numoutside)
    // for (int i = 0; i < npoints; i++) {
    //     for (int j = 0; j < npoints; j++) {
    //         // Define c for each point
    //         c_real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
    //         c_imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
            
    //         // Initialize z for each point
    //         double z_real = c_real;
    //         double z_imag = c_imag;

    //         // Perform iterations for each point
    //         for (int iter = 0; iter < maxiter; iter++) {
    //             ztemp_real = (z_real * z_real) - (z_imag * z_imag) + c_real;
    //             ztemp_imag = z_real * z_imag * 2 + c_imag;
    //             z_real = ztemp_real;
    //             z_imag = ztemp_imag;
                
    //             // Check if z is outside the boundary
    //             if ((z_real * z_real + z_imag * z_imag) > 4.0e0) {
    //                 numoutside++;
    //                 break;
    //             }
    //         }
    //     }
    // }

    // double ztemp_real, ztemp_imag;
    // double c_real, c_imag;

    // // Parallelize the nested loops and use reduction for numoutside
    // #pragma omp parallel for collapse(2) private(c_real, c_imag, ztemp_real, ztemp_imag) reduction(+:numoutside)
    // for (int i = 0; i < npoints; i++) {
    //     for (int j = 0; j < npoints; j++) {
    //         // Define c for each point
    //         c_real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
    //         c_imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
            
    //         // Initialize z for each point
    //         double z_real = c_real;
    //         double z_imag = c_imag;

    //         // Perform iterations for each point
    //         for (int iter = 0; iter < maxiter; iter++) {
    //             ztemp_real = (z_real * z_real) - (z_imag * z_imag) + c_real;
    //             ztemp_imag = z_real * z_imag * 2 + c_imag;
    //             z_real = ztemp_real;
    //             z_imag = ztemp_imag;
                
    //             // Check if z is outside the boundary
    //             if ((z_real * z_real + z_imag * z_imag) > 4.0e0) {
    //                 numoutside++;
    //                 break;
    //             }
    //         }
    //     }
    // }

    // Correct code, but slow
    // double ztemp_real, ztemp_imag;
    // double c_real, c_imag;

    // // Parallelize the outer loop
    // #pragma omp parallel for reduction(+:numoutside)
    // for (int i = 0; i < npoints; i++) {
    //     // Parallelize the inner loop
    //     #pragma omp parallel for reduction(+:numoutside)
    //     for (int j = 0; j < npoints; j++) {
    //         // Define c for each point
    //         c_real = -2.0 + 2.5 * (double)(i) / (double)(npoints) + 1.0e-7;
    //         c_imag = 1.125 * (double)(j) / (double)(npoints) + 1.0e-7;
            
    //         // Initialize z for each point
    //         double z_real = c_real;
    //         double z_imag = c_imag;

    //         // Perform iterations for each point
    //         for (int iter = 0; iter < maxiter; iter++) {
    //             ztemp_real = (z_real * z_real) - (z_imag * z_imag) + c_real;
    //             ztemp_imag = z_real * z_imag * 2 + c_imag;
    //             z_real = ztemp_real;
    //             z_imag = ztemp_imag;
                
    //             // Check if z is outside the boundary
    //             if ((z_real * z_real + z_imag * z_imag) > 4.0e0) {
    //                 numoutside++;
    //                 break;
    //             }
    //         }
    //     }
    // }

    #pragma omp parallel for default(none) shared(npoints, maxiter)  private(i, j, c, z, iter, ztemp) schedule(dynamic) reduction(+:numoutside)
    for (int i = 0; i < npoints; i++) {
        for (int j = 0; j < npoints; j++) {
            //complex c;
            c.real = -2.0 + 2.5 * (double)i / (double)npoints + 1.0e-7;
            c.imag = 1.125 * (double)j / (double)npoints + 1.0e-7;

            //complex z = c;
            z = c;
            for (int iter = 0; iter < maxiter; iter++) {
                double ztemp = (z.real * z.real) - (z.imag * z.imag) + c.real;
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
