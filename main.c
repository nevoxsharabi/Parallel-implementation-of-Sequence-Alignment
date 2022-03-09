#include "mpi.h"
#include "omp.h"
#include "mpiUtil.h"

float findSimilarityWeight(char a, char b, float* weights);
int isIdentical(char a, char b);
int is_conservative(char a, char b);
int is_semi_conservative(char a, char b);

int main(int argc, char* argv[])
{
	int rank, num_procs, numOfSeq2s = 0;
	char seq1[SEQ1_LEN], seq2[SEQ2_LEN];
	float weights[WEIGHTS_NUM];
	Bundle bundle;	
	
	MPI_Datatype resultType;
	MPI_Datatype bundleType;
	MPI_Status status;

	InitMPI(&argc, &argv, &rank, &num_procs, &resultType, &bundleType);

	if (num_procs != NUM_PROCS)
	{
		printf("please use this program with %d processes\n", NUM_PROCS);
		MPI_Abort(MPI_COMM_WORLD, __LINE__);
	}
	
	int start_time = omp_get_wtime();
	
	if (rank == MASTER) 
	{
		char** seqs2 = readFromFile(READ_FILE_NAME, weights, seq1, &numOfSeq2s); // get input data
		printf("number of seq2 is: %d\n", numOfSeq2s);
		// Creating Conservative Matrix
		int size_of_cons_matrix = 26 * 26;
		float conservative_matrix[size_of_cons_matrix];
		for(int i = 0; i < 26; i++){
    		   for(int j = 0; j < 26; j++){
     		      conservative_matrix[i * 26 + j] = findSimilarityWeight(i + 'A', j + 'A', &weights[0]);
     		      //printf("%1.2f ", conservative_matrix[i * 26 + j]);
    		   }
    		   //printf("\n");
   		}
		
		for(int i = 0; i < numOfSeq2s; i++){
			printf("working on SEQ2_%d\n", i);
			calcBestScoreCUDA(seq1, seqs2[i], &conservative_matrix[0], &weights[0]);
		}
		
	}
	else 
	{

	}
	
	printf("process %d finished work after - %.2f seconds\n", rank, omp_get_wtime() - start_time);
		
	MPI_Finalize();
	return 0;
}

float findSimilarityWeight(char a, char b, float* weights){
    if (isIdentical(a, b))
        return weights[0];
    else if(is_conservative(a, b))
        return -weights[1];
    else if(is_semi_conservative(a, b))
        return -weights[2];
    else
        return -weights[3];
}

int isIdentical(char a, char b){
    if(a == b)
        return 1;
    return 0;
}

int is_conservative(char a, char b){
    const char* conservative_groups[9] = {"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"};
	for(int i = 0; i < 9; i++){
		if(strchr((char*)conservative_groups[i], a) && strchr((char*)conservative_groups[i], b))
			return 1;
	}
	return 0;
}

int is_semi_conservative(char a, char b){
    const char* semi_conservative_groups[11] = {"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"};
	for(int i = 0; i < 11; i++){
		if(strchr((char*)semi_conservative_groups[i], a) && strchr((char*)semi_conservative_groups[i], b))
			return 1;
	}
	return 0;
}
