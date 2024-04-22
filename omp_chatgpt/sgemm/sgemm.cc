#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <malloc.h>
#include <vector>
#include <iostream>
#include<fstream>

#include <omp.h>

bool readColMajorMatrixFile(const char *fn, int &nr_row, int &nr_col, std::vector<float>&v)
{
  std::cerr << "Opening file:"<< fn << std::endl;
  std::fstream f(fn, std::fstream::in);
  if ( !f.good() ) {
    return false;
  }

  // Read # of rows and cols
  f >> nr_row;
  f >> nr_col;

  float data;
  std::cerr << "Matrix dimension: "<<nr_row<<"x"<<nr_col<<std::endl;
  while (f.good() ) {
    f >> data;
    v.push_back(data);
  }
  v.pop_back(); // remove the duplicated last element
  return true;
}

bool writeColMajorMatrixFile(const char *fn, int nr_row, int nr_col, std::vector<float>&v)
{
  std::cerr << "Opening file:"<< fn << " for write." << std::endl;
  std::fstream f(fn, std::fstream::out);
  if ( !f.good() ) {
    return false;
  }

  // Read # of rows and cols
  f << nr_row << " "<<nr_col<<" ";

  float data;
  std::cerr << "Matrix dimension: "<<nr_row<<"x"<<nr_col<<std::endl;
  for (int i = 0; i < v.size(); ++i) {
    f << v[i] << ' ';
  }
  f << "\n";
  return true;

}

/* 
 * Base C implementation of MM
 */
/*
Apologies for the oversight. Let's correct the code. The issue arises because we're trying to perform arithmetic operations inside the #pragma omp atomic directive, which isn't allowed. We'll adjust the code to first calculate the result outside of the atomic region and then perform the update.

void basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if ((transa != 'N') && (transa != 'n'))
    {
        std::cerr << "unsupported value of 'transa'" << std::endl;
        return;
    }

    if ((transb != 'T') && (transb != 't'))
    {
        std::cerr << "unsupported value of 'transb'" << std::endl;
        return;
    }

    #pragma omp parallel for collapse(2) // collapse the nested loops for better parallelization
    for (int mm = 0; mm < m; ++mm)
    {
        for (int nn = 0; nn < n; ++nn)
        {
            float c = 0.0f;
            for (int i = 0; i < k; ++i)
            {
                float a = A[mm + i * lda];
                float b = B[nn + i * ldb];
                c += a * b;
            }
            #pragma omp atomic // ensure atomic update for concurrent writes to C
            C[mm + nn * ldc] = C[mm + nn * ldc] * beta; //+ alpha * c;
        }
    }
}
*/
/*
I see, my apologies for the misunderstanding. The #pragma omp atomic directive is not appropriate in this context because it's designed for performing simple atomic operations like increments or assignments, but not for more complex expressions.
void basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if ((transa != 'N') && (transa != 'n'))
    {
        std::cerr << "unsupported value of 'transa'" << std::endl;
        return;
    }

    if ((transb != 'T') && (transb != 't'))
    {
        std::cerr << "unsupported value of 'transb'" << std::endl;
        return;
    }

    #pragma omp parallel for collapse(2)
    for (int mm = 0; mm < m; ++mm)
    {
        for (int nn = 0; nn < n; ++nn)
        {
            float c = 0.0f;
            for (int i = 0; i < k; ++i)
            {
                float a = A[mm + i * lda];
                float b = B[nn + i * ldb];
                c += a * b;
            }
            float result = C[mm + nn * ldc] * beta + alpha * c;
            #pragma omp atomic
            C[mm + nn * ldc] = result;
        }
    }
}
*/

void basicSgemm(char transa, char transb, int m, int n, int k, float alpha, const float *A, int lda, const float *B, int ldb, float beta, float *C, int ldc)
{
    if ((transa != 'N') && (transa != 'n'))
    {
        std::cerr << "unsupported value of 'transa'" << std::endl;
        return;
    }

    if ((transb != 'T') && (transb != 't'))
    {
        std::cerr << "unsupported value of 'transb'" << std::endl;
        return;
    }

    #pragma omp parallel for collapse(2)
    for (int mm = 0; mm < m; ++mm)
    {
        for (int nn = 0; nn < n; ++nn)
        {
            float c = 0.0f;
            for (int i = 0; i < k; ++i)
            {
                float a = A[mm + i * lda];
                float b = B[nn + i * ldb];
                c += a * b;
            }
            C[mm + nn * ldc] = C[mm + nn * ldc] * beta + alpha * c;
        }
    }
}

int main (int argc, char *argv[]) {

  int matArow, matAcol;
  int matBrow, matBcol;
  std::vector<float> matA, matBT;
  double timer_start, timer_end, elapsed_time;

  if (argc != 4)  
  {
      fprintf(stderr, "Expecting three input filenames\n");
      exit(-1);
  }
 
  /* Read in data */
  // load A
  readColMajorMatrixFile(argv[1], matArow, matAcol, matA);

  // load B^T
  readColMajorMatrixFile(argv[2], matBcol, matBrow, matBT);

  // allocate space for C
  std::vector<float> matC(matArow*matBcol);

  // Use standard sgemm interface
  timer_start = omp_get_wtime();
  basicSgemm('N', 'T', matArow, matBcol, matAcol, 1.0f, &matA.front(), matArow, &matBT.front(), matBcol, 0.0f, &matC.front(), matArow);
  timer_end = omp_get_wtime();
  writeColMajorMatrixFile(argv[3], matArow, matBcol, matC); 

  elapsed_time = timer_end - timer_start;
  printf("Elapsed time: %.6f\n", elapsed_time);

  return 0;
}
