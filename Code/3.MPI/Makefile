#TODO: add the following implementations: multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

all: hooke

hooke: multistart_hooke_mpi.c
	mpicc -Wall -O2 -o multistart_hooke_mpi multistart_hooke_mpi.c -lm

mac: multistart_hooke_mpi.c
	mpicc -Wall -O2 -o multistart_hooke_mpi multistart_hooke_mpi.c -lm


clean:
	rm -f multistart_hooke_mpi
