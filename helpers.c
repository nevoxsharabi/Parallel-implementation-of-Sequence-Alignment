#include "header.h"

// ========================

void assignOffsets(int* start, int* end, int maxOffset, int rank) // divide each process portion of offsets
{
	int portion = maxOffset / NUM_PROCS;
	*start = rank * portion;

	// if offset is odd give remainder to last process
	if (maxOffset % NUM_PROCS != 0 && rank == NUM_PROCS - 1)
		portion += maxOffset % NUM_PROCS;

	*end = *start + portion;
}

// ========================

char* MS(char* seq, int n, int k) // remove chars in n & k index
{
	int i = n, j = k-1; 
	
	char* mutant = (char*)malloc(strlen(seq) * sizeof(char));
	
	strcpy(mutant,seq);
	
	memmove(&mutant[i], &mutant[i + 1], strlen(mutant) - i);
	memmove(&mutant[j], &mutant[j + 1], strlen(mutant) - j);
	
	return mutant;
}

// ========================

int** createNKs(int num_mutants) // create all n k for mutants
{
	int** nkArr = (int**) malloc(num_mutants * sizeof(int*));

	int n = 0, k = 1;
	for (int i = 0; i < num_mutants; i++) 
	{
		nkArr[i] = (int*) malloc(2*sizeof(int));
		nkArr[i][0] = n;
		nkArr[i][1] = k;

		n++;
		if (n == k)
		{
			n = 0;
			k++;
		} 
	}
	
	return nkArr;
}

// ========================
