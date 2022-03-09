#include "mpiUtil.h"
#include <stddef.h>

#define RESULT_ATTRS_NUM 4
#define BUNDLE_ATTRS_NUM 5

// ========================

// initialize mpi communication
void InitMPI(int* argc, char** argv[], int* rank, int* num_procs, MPI_Datatype* resultType, MPI_Datatype* bundleType)
{
	// Init MPI
	MPI_Init(argc, argv);
	MPI_Comm_rank(MPI_COMM_WORLD, rank); // get process rank
	MPI_Comm_size(MPI_COMM_WORLD, num_procs); // get number of processors

	createResultType(resultType); //create resultType
	createBundleType(bundleType); //create bundleType
}

// ========================

// create data type for mpi communication
void createResultType(MPI_Datatype* resultType)
{
	int blockLen[RESULT_ATTRS_NUM] = { 1, 1, 1, 1};
	MPI_Aint disp[RESULT_ATTRS_NUM];
	MPI_Datatype types[RESULT_ATTRS_NUM] = { MPI_FLOAT, MPI_INT, MPI_INT, MPI_INT };

	disp[0] = offsetof(Result, score);
	disp[1] = offsetof(Result, offset);
	disp[2] = offsetof(Result, n);
	disp[3] = offsetof(Result, k);

	MPI_Type_create_struct(RESULT_ATTRS_NUM, blockLen, disp, types, resultType);
	MPI_Type_commit(resultType);
}

// ========================

// create data type for data bundle
void createBundleType(MPI_Datatype* bundleType)
{
	int blockLen[BUNDLE_ATTRS_NUM] = { SEQ1_LEN, 1, SEQ2_LEN, 1, WEIGHTS_NUM};
	MPI_Aint disp[BUNDLE_ATTRS_NUM];
	MPI_Datatype types[BUNDLE_ATTRS_NUM] = { MPI_CHAR, MPI_INT, MPI_CHAR, MPI_INT, MPI_FLOAT };

	disp[0] = offsetof(Bundle, seq1);
	disp[1] = offsetof(Bundle, seq1_len);
	disp[2] = offsetof(Bundle, seq2);
	disp[3] = offsetof(Bundle, seq2_len);
	disp[4] = offsetof(Bundle, weights);
	
	MPI_Type_create_struct(BUNDLE_ATTRS_NUM, blockLen, disp, types, bundleType);
	MPI_Type_commit(bundleType);
}

// ========================

void populateBundleData(char* seq1, char* seq2, float weights[], Bundle* bundle)
{
	int len = strlen(seq1);
	for (int i = 0; i < len; i++)
		bundle->seq1[i] = seq1[i];
	bundle->seq1_len = len;
	
	len = strlen(seq2);
	for (int i = 0; i < len; i++)
		bundle->seq2[i] = seq2[i];
	bundle->seq2_len = len;
	
	for (int i = 0; i < WEIGHTS_NUM; i++)
		bundle->weights[i] = weights[i];	
}
