/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i;
    int j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(i, j, prime) reduction(+:total)
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
        total += prime;
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
}*/

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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

        // Parallelize the outer loop
        // error: invalid increment expression
        #pragma omp parallel for reduction(+:primes)
        for (int i = n_lo; i <= n; i *= n_factor)
        {
            primes = prime_number(i);
            printf("  %8d  %8d  %14f\n", i, primes, omp_get_wtime() - ctime);
        }

        n *= n_factor;
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
*/

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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

        // Parallelize the outer loop
        #pragma omp parallel for reduction(+:primes)
        for (int i = n_lo; i <= n; i++)
        {
            primes = prime_number(i);
            printf("  %8d  %8d  %14f\n", i, primes, omp_get_wtime() - ctime);
        }

        n *= n_factor;
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
}*/

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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

        // Parallelize the outer loop
        #pragma omp parallel for reduction(+:primes) firstprivate(n)
        for (int i = n_lo; i <= n_hi; i += n_factor)
        {
            primes = prime_number(n);
            printf("  %8d  %8d  %14f\n", n, primes, omp_get_wtime() - ctime);
            n *= n_factor; // Calculate the next value of n
            if (n > n_hi) break; // Avoid exceeding n_hi
        }
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
}*/

/* reduntant iteratios */
/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
    double ctime;
    int stop = 0; // Variable to control loop termination

    printf("\n");
    printf("  Call PRIME_NUMBER to count the primes from 1 to N.\n");
    printf("\n");
    printf("         N        Pi          Time\n");
    printf("\n");

    n = n_lo;

    while (!stop)
    {
        ctime = omp_get_wtime();

        // Parallelize the outer loop
        #pragma omp parallel for reduction(+:primes) firstprivate(n)
        for (int i = n_lo; i <= n_hi; i += n_factor)
        {
            if (i <= n_hi)
            {
                primes = prime_number(n);
                printf("  %8d  %8d  %14f\n", n, primes, omp_get_wtime() - ctime);
                n *= n_factor; // Calculate the next value of n
            }
            else
            {
                // If the next iteration exceeds n_hi, set the termination flag
                #pragma omp atomic write
                stop = 1;
            }
        }
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
}*/

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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

        // Parallelize the outer loop
        #pragma omp parallel for reduction(+:primes) firstprivate(n)
        for (int i = 0; i < n_factor && n <= n_hi; i++)
        {
            primes = prime_number(n);
            printf("  %8d  %8d  %14f\n", n, primes, omp_get_wtime() - ctime);
            n += n_lo; // Increment n by n_lo
        }
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
*/

/*#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        for (j = 2; j * j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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
        printf("  %8d  %8d  %14f\n", n, primes, omp_get_wtime() - ctime);
        n *= n_factor;
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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
}*/

/* Missing test function */
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total) schedule(dynamic, 1000)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
//        for (j = 2; j * j <= i; j++)
        for (j = 2; j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

int main(int argc, char *argv[])
{
    int n = 1000000;
    int num_primes;

    if (argc > 1)
        n = atoi(argv[1]);

    double start_time = omp_get_wtime();
    num_primes = prime_number(n);
    double end_time = omp_get_wtime();

    printf("Number of primes up to %d: %d\n", n, num_primes);
    printf("Execution time: %f seconds\n", end_time - start_time);

    return 0;
}*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int prime_number(int n)
{
    int i, j;
    int prime;
    int total = 0;

    #pragma omp parallel for private(j, prime) reduction(+:total) schedule(dynamic, 1000)
    for (i = 2; i <= n; i++)
    {
        prime = 1;
        // Optimized loop: only check up to square root of i
        // for (j = 2; j * j <= i; j++)
        for (j = 2; j <= i; j++)
        {
            if ((i % j) == 0)
            {
                prime = 0;
                break;
            }
        }
        total += prime;
    }
    return total;
}

void test(int n_lo, int n_hi, int n_factor)
{
    int n, primes;
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
        printf("  %8d  %8d  %14f\n", n, primes, omp_get_wtime() - ctime);
        n *= n_factor;
    }
}

int main(int argc, char *argv[])
{
    int n_factor, n_hi, n_lo;

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






