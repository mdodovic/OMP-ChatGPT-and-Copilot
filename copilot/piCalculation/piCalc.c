#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void Usage(char* prog_name);

int main(int argc, char* argv[]) {
   long long n, i;
   double factor;
   double sum = 0.0;
   double timer_start, timer_end, elapsed_time;


   if (argc != 2) Usage(argv[0]);
   n = strtoll(argv[1], NULL, 10);
   if (n < 1) Usage(argv[0]);

   printf("Before for loop, factor = %f.\n", factor);
   timer_start = omp_get_wtime();

# pragma omp parallel for reduction(+:sum)
   for (i = 0; i < n; i++) {
      factor = (i % 2 == 0) ? 1.0 : -1.0; 
      sum += factor/(2*i+1);
   }
   timer_end = omp_get_wtime();
   printf("After for loop, factor = %f.\n", factor);
   elapsed_time = timer_end - timer_start;
   
   sum = 4.0*sum;
   printf("With n = %lld terms\n", n);
   printf("   Our estimate of pi = %.14f\n", sum);
   printf("   Ref estimate of pi = %.14f\n", 4.0*atan(1.0));
   printf("Elapsed time: %.6f\n", elapsed_time);
   return 0;
}

void Usage(char* prog_name) {
   fprintf(stderr, "usage: %s <thread_count> <n>\n", prog_name);
   fprintf(stderr, "   n is the number of terms and should be >= 1\n");
   exit(0);
}
