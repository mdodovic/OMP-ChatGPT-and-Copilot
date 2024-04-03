# C compiler
CC = gcc
CC_FLAGS = -std=gnu99 -fopenmp -O2

hotspot_for: 
	$(CC) $(CC_FLAGS) hotspot_for.c -o hotspot_for 

hotspot_task: 
	$(CC) $(CC_FLAGS) hotspot_task.c -o hotspot_task 

clean:
	rm -f hotspot_for hotspot_task
