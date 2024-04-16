#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MINREAL -1024.0
#define MAXREAL 1024.0

#include <sys/time.h>

double gettime(void) {
   struct timeval tv;
   gettimeofday(&tv, NULL);
   return tv.tv_sec + 1e-6 * tv.tv_usec;
}

void matFillSimple(long N, double *mat, double val) {
   long i, j;

   for(i = 0; i < N; i ++)
      for(j = 0; j < N; j ++)
         mat[i * N + j] = val;
}

void matFillRand(long N, double *mat, double val) {
   long i, j;

   for(i = 0; i < N; i ++)
      for(j = 0; j < N; j ++)
         mat[i * N + j] = (rand() / (double) RAND_MAX)*(MAXREAL - MINREAL) + MINREAL;
}

void matFillTriangle(long N, double *mat, double val) {
   long i, j;

   for(i = 0; i < N; i ++)
      for(j = 0; j < N; j ++) { 
         if (i < j) mat[i * N + j] = (rand() / (double) RAND_MAX)*(MAXREAL - MINREAL) + MINREAL;
         else  mat[i * N + j] = 0.0;
      }
}

void matMul(long N, double *a, double *b, double *c) {
   long i, j, k;

   for (i = 0; i < N; i ++)
      for (j = 0; j < N; j ++)
         for (k = 0; k < N; k ++)
            c[i * N + j] += a[i * N + k] * b[k * N + j];
}

int main(int argc, char **argv) {
   long N;
   double *A, *B, *C, t;

   srand(time(NULL));

   N = atoi(argv[1]);
   A = (double *) malloc(N * N * sizeof(double));
   B = (double *) malloc(N * N * sizeof(double));
   C = (double *) malloc(N * N * sizeof(double));
   matFillSimple(N, A, 1.0);
   matFillSimple(N, B, 2.0);
   matFillSimple(N, C, 0.0);

   t = gettime();
   matMul(N, A, B, C);
   t = gettime() - t;

   fprintf(stdout, "%ld\t%lf\n", N, t);
   fflush(stdout);

   free(A);
   free(B);
   free(C);

   return EXIT_SUCCESS;
}
