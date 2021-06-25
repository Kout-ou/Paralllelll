#TODO: add the following implementations: multistart_hooke_omp multistart_hooke_omp_tasks multistart_hooke_mpi multistart_hooke_mpi_omp

all: hooke

hooke: multistart_hooke_hybrid.c
	 mpicc -Wall -O3 -o multistart_hooke_hybrid multistart_hooke_hybrid.c -lm -fopenmp

mac: multistart_hooke_hybrid.c
	mpicc -Wall -O3 -o multistart_hooke_hybrid multistart_hooke_hybrid.c -lm -openmp

clean:
	rm -f multistart_hooke_mpi_omp
