# C compiler
CC = gcc
CC_FLAGS = -fopenmp -O3

feyman: 
	$(CC) $(CC_FLAGS) feyman.c -o feyman -lm

feyman_omp: 
	$(CC) $(CC_FLAGS) feyman_omp.c -o feyman_omp -lm

clean:
	rm -f feyman feyman_omp
