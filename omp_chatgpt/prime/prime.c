#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <omp.h>

int prime_number1(int n) {
    int total = 0;

    #pragma omp parallel for reduction(+:total)
    for (int i = 2; i <= n; i++) {
        int prime = 1;
        for (int j = 2; j < i; j++) {
            if ((i % j) == 0) {
                prime = 0;
                break;
            }
        }
        total += prime;
    }

    return total;
}

// ChatGPT offered optimized version with changed algorithm
// In this case loop condition has been changed to j * j <= i
int prime_number(int n) {
    int total = 0;

    #pragma omp parallel for reduction(+:total) schedule(dynamic)
    for (int i = 2; i <= n; i++) {
        int prime = 1;
        //for (int j = 2; j * j <= i; j++) {
        for (int j = 2; j <= i; j++) {
            if (i % j == 0) {
                prime = 0;
                break;
            }
        }
        total += prime;
    }

    return total;
}

int prime_number_s(int n)
{
  int i;
  int j;
  int prime;
  int total;

  total = 0;

  for (i = 2; i <= n; i++)
  {
    prime = 1;
    for (j = 2; j < i; j++)
    {
      if ((i % j) == 0)
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }
  return total;
}

void test(int n_lo, int n_hi, int n_factor);

int main(int argc, char *argv[])
{
  int n_factor;
  int n_hi;
  int n_lo;

  printf("PRIME TEST\n");

  if (argc != 4)
  {
    n_lo = 1;
    n_hi = 131072;
    n_factor = 2;
  }
  else
  {
    n_lo = atoi(argv[1]);
    n_hi = atoi(argv[2]);
    n_factor = atoi(argv[3]);
  }

  test(n_lo, n_hi, n_factor);

  printf("\n");
  printf("PRIME_TEST\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  
  return 0;
}

void test(int n_lo, int n_hi, int n_factor)
{
  int i;
  int n;
  int primes;
  double ctime;

  printf("\n");
  printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
  printf("\n");
  printf("         N        Pi          Time\n");
  printf("\n");

  n = n_lo;

  while (n <= n_hi)
  {
    ctime = omp_get_wtime();

    primes = prime_number(n);

    ctime = omp_get_wtime() - ctime;

    printf("  %8d  %8d  %14f\n", n, primes, ctime);
    n = n * n_factor;
  }

  return;
}
