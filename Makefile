build:
	mpicxx -fopenmp -c main.c -o main.o
	mpicxx -fopenmp -c fileUtil.c -o fileUtil.o
	mpicxx -fopenmp -c mpiUtil.c -o mpiUtil.o
	mpicxx -fopenmp -c helpers.c -o helpers.o
	mpicxx -fopenmp -c cFunctions.c -o cFunctions.o
	nvcc -I./common/inc -c cudaFunctions.cu -o cudaFunctions.o
	mpicxx -fopenmp -o program main.o fileUtil.o mpiUtil.o helpers.o cFunctions.o cudaFunctions.o /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt
clean:
	rm -f *.o ./program
	
run:
	mpiexec -np 2 ./program
