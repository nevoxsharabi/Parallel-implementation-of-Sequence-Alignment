#include "header.h"
#include "omp.h"
#include <float.h>

void printBundle(Bundle bundle)
{
	printf("%s\n", bundle.seq1);
	printf("%s\n", bundle.seq2);
	printf("%d\n", bundle.seq1_len);
	printf("%d\n", bundle.seq2_len);
	for (int i = 0; i < WEIGHTS_NUM; i++)
		printf("%.1f\n", bundle.weights[i]);
}

// ========================

Result calcBestScore(Bundle bundle, int rank)
{
	int len1 = bundle.seq1_len;
	int len2 = bundle.seq2_len;
	int maxOffset = len1 - (len2-2);
	int num_mutants = len2 * (len2 - 1) / 2;
	
	Result res;
	res.score = -FLT_MAX;
	
	#pragma omp parallel for
	for (int n = 0; n < len2; n++) {
		//for (int k = n+1; k < len2; k++) 
	}
	
	return res;
}

// ========================

Result calcOffsetBestResult(char* seq1, char* seq2, float weights[], int num_mutants, int offset)
{
	int len = strlen(seq2);
	int n = 0, k = 1;

	Result bestRes;
	bestRes.score = -FLT_MAX;
	
	int** nkArr = createNKs(num_mutants); // create n k for all mutants
	
	#pragma omp parallel for num_threads(4) schedule(static, num_mutants/4)
	for (int i = 0; i < num_mutants; i++)
	{
		char* mutant = MS(seq2, nkArr[i][0], nkArr[i][1]); // create mutant sequence
		float* similarity = (float*)malloc((len - 2) * sizeof(float));
		
		calcMutantSimilarityCUDA(similarity, seq1 + offset, mutant, weights); // compare seq1 and mutant and assign weights
		
		Result res = calcMutantResultOMP(similarity, len - 2, offset, nkArr[i][0], nkArr[i][1]); // sum weights and find mutant score

		bestRes = res.score > bestRes.score ? res : bestRes; // swap if found better score
		
		free(similarity);
	}

	return bestRes;
}

// ========================

Result calcMutantResultOMP(float* similarity, int mutant_size, int offset, int n, int k)
{
	Result res;
	float score = 0;
	
	#pragma omp parallel for reduction(+:score)
	for(int i = 0; i < mutant_size; i++) // sum mutant weights
		score += similarity[i]; 
	
	//create result struct for mutant
	res.score = score;
	res.offset = offset;
	res.n = n+1;
	res.k = k+1;

	return res;
}

// ========================

void calcBestScoreOmp(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants){
	float sscore[4] = {-100000, -100000, -100000, -100000};
	int ooffset[4];
	int mmutant[4];
	#pragma omp parallel for
	for(int i = 0; i < num_mutants; i++){
		int tid = omp_get_thread_num();
		if(mutantsBestScores[i] > sscore[tid]){
			sscore[tid] = mutantsBestScores[i];
			ooffset[tid] = mutantsBestOffsets[i];
			mmutant[tid] = i;
			
		}
	}
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[0], ooffset[0], mmutant[0]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[1], ooffset[1], mmutant[1]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[2], ooffset[2], mmutant[2]);
	//printf("score = %1.2f, offset = %d, mutant_num = %d\n", sscore[3], ooffset[3], mmutant[3]);
	float best_score = -100000;
	int best_offset;
	int mutant_num;
	for(int i = 0; i < 4; i++){
		if(sscore[i] > best_score){
			best_score = sscore[i];
			best_offset = ooffset[i];
			mutant_num = mmutant[i];
		}
	}
	printf("best score = %1.2f, best offset = %d, mutant num = %d\n", best_score, best_offset, mutant_num);
}

