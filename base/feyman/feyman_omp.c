#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>

int i4_ceiling(double x) {
    int value = (int)x;
    if (value < x)
        value = value + 1;
    return value;
}

int i4_min(int i1, int i2) {
    return (i1 < i2) ? i1 : i2;
}

double potential(double a, double b, double c, double x, double y, double z) {
    return 2.0 * (pow(x / a / a, 2) + pow(y / b / b, 2) + pow(z / c / c, 2)) + 1.0 / a / a + 1.0 / b / b + 1.0 / c / c;
}

double r8_uniform_01(int *seed) {
    int k;
    double r;

    k = *seed / 127773;

    *seed = 16807 * (*seed - k * 127773) - k * 2836;

    if (*seed < 0) {
        *seed = *seed + 2147483647;
    }
    r = (double)(*seed) * 4.656612875E-10;

    return r;
}

void timestamp(void) {
    #define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    printf("%s\n", time_buffer);

    #undef TIME_SIZE
}

int main(int argc, char **argv) {
    double a = 3.0;
    double b = 2.0;
    double c = 1.0;
    double err = 0.0;
    double h = 0.001;
    int n_inside = 0;
    int ni, nj, nk;
    int N = atoi(argv[1]);

    timestamp();

    printf("A = %f\n", a);
    printf("B = %f\n", b);
    printf("C = %f\n", c);
    printf("N = %d\n", N);
    printf("H = %6.4f\n", h);

    ni = 6;
    nj = 1 + i4_ceiling(b / a) * (ni - 1);
    nk = 1 + i4_ceiling(c / a) * (ni - 1);

    err = 0.0;

    #pragma omp parallel for reduction(+:n_inside, err)
    for (int i = 1; i <= ni; i++) {
        double x = ((double)(ni - i) * (-a) + (double)(i - 1) * a) / (double)(ni - 1);
        for (int j = 1; j <= nj; j++) {
            double y = ((double)(nj - j) * (-b) + (double)(j - 1) * b) / (double)(nj - 1);
            for (int k = 1; k <= nk; k++) {
                double z = ((double)(nk - k) * (-c) + (double)(k - 1) * c) / (double)(nk - 1);
                double chk = pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2);
                if (1.0 < chk) {
                    continue;
                }
                n_inside++;
                double w_exact = exp(pow(x / a, 2) + pow(y / b, 2) + pow(z / c, 2) - 1.0);
                double wt = 0.0;
                for (int trial = 0; trial < N; trial++) {
                    double x1 = x, x2 = y, x3 = z;
                    double w = 1.0, chk = 0.0;
                    while (chk < 1.0) {
                        double dx = ((double)rand() / RAND_MAX - 0.5) * sqrt(3.0 * h);
                        double dy = ((double)rand() / RAND_MAX - 0.5) * sqrt(3.0 * h);
                        double dz = ((double)rand() / RAND_MAX - 0.5) * sqrt(3.0 * h);
                        double vs = potential(a, b, c, x1, x2, x3);
                        x1 += dx; x2 += dy; x3 += dz;
                        double vh = potential(a, b, c, x1, x2, x3);
                        double we = (1.0 - h * vs) * w;
                        w -= 0.5 * h * (vh * we + vs * w);
                        chk = pow(x1 / a, 2) + pow(x2 / b, 2) + pow(x3 / c, 2);
                    }
                    wt += w;
                }
                wt /= (double)N;
                err += pow(w_exact - wt, 2);
            }
        }
    }

    err = sqrt(err / (double)n_inside);

    printf("\n\nRMS absolute error in solution = %e\n", err);
    timestamp();

    return 0;
}
