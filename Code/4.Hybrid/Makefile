#TODO: add the following implementations: multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

all: multistart_hooke_seq

hooke: multistart_hooke_seq.c
	 gcc -Wall -O3 -o multistart_hooke_mpi_omp multistart_hooke_mpi_omp.c -lm

clean:
	rm -f multistart_hooke_mpi_omp
