#TODO: add the following implementations: multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

all: hooke

hooke: multistart_hooke_omp_tasks.c
	 gcc -Wall -O3 -o multistart_hooke_omp_tasks multistart_hooke_omp_tasks.c -lm

mac: multistart_hooke_omp_tasks.c
	gcc-11 -Wall -O3 -o multistart_hooke_omp_tasks multistart_hooke_omp_tasks.c -lm

clean:
	rm -f multistart_hooke_omp_tasks
