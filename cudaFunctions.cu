#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "header.h"

// ======================== macros for allocating device memory and error check

#define CUDA_ERR_CHECK(err,msg) (\
		{if (err != cudaSuccess) { \
			fprintf(stderr, msg " - %s\n", cudaGetErrorString(err)); \
			exit(EXIT_FAILURE); \
		} \
	})

#define CUDA_MEM_INIT(dest, src, size, type) {\
	cudaError_t err = cudaSuccess;\
	size_t  arrSize = size * sizeof(type);\
	err = cudaMalloc((void**)&dest, arrSize);\
	CUDA_ERR_CHECK(err, "Failed to allocate device memory");\
	err = cudaMemcpy(dest, src, arrSize, cudaMemcpyHostToDevice);\
	CUDA_ERR_CHECK(err, "Failed to copy data from host to device"); }\

	
// ======================== groups
void getNK_CPU(int mutant_num, int seq2_len, int* n, int* k);
void getNK_CPU(int mutant_num, int seq2_len, int* n, int* k)
{
	int i;
	int num_of_mutants_in_row = seq2_len;

	for(i = 1; i < seq2_len; i++){
		if(mutant_num - (num_of_mutants_in_row - 1) > 0){
		    mutant_num -= (num_of_mutants_in_row - 1);
		    num_of_mutants_in_row--;
		}else{
		    break;
		}
	}
	
	*n = i;	
	*k = i + mutant_num;
}

__device__ const char* conservatives[CONSERVATIVES_LEN] = {
	"NDEQ", "NEQK", "STA", "MILV", "QHRK", "NHQK", "FYW", "HY", "MILF"
};

__device__ const char* semi_conservatives[SEMI_CONSERVATIVES_LEN] = {
	"SAG", "ATV", "CSA", "SGND", "STPA", "STNK", "NEQHRK", "NDEQHK", "SNDEQK", "HFY", "FVLIM"
};

// ======================== // check if char in group (same as strchr)

__device__ int isCharExist(const char* s, int c)
{
	do {
	    if (*s == c)
		return 1;
	} while (*s++);
	  
  	return 0;
}

// ========================

__device__ int checkInGroup(char first, char second, const char** group, int len) // checks if both chars in the group
{
	for (int i = 0; i < len; i++)
	{
		if (isCharExist(group[i], first) != NULL
			&& isCharExist(group[i], second) != NULL)
			return 1;
	}
	return 0;
}

// ========================

__device__ float compareChars(char first, char second, float weights[]) // assign suitable weight for chars 
{
	if (first == second)
		return weights[0];
	else if (checkInGroup(first, second, conservatives, CONSERVATIVES_LEN))
		return -weights[1];
	else if (checkInGroup(first, second, semi_conservatives, SEMI_CONSERVATIVES_LEN))
		return -weights[2];
		
	return -weights[3];
}

// ========================

__global__ void calcSimilarityKernel(float* d_similarity, char* d_seq1, char* d_mutant, float* d_weights, int mutant_size) // traverses char by char and checks similarity
{
	//int i = blockIdx.x * blockDim.x + threadIdx.x;
	
	//if(i < mutant_size)
		//d_similarity[i] = compareChars(d_seq1[i], d_mutant[i], d_weights);
}

// ========================

void calcMutantSimilarityCUDA(float* similarity, char* seq1, char* mutant, float weights[]) // entry point to CUDA , allocates memory and calls kernel func
{
	int mutant_size = strlen(mutant);
	
	char* d_seq1 = NULL; // allocate seq1 memory
	CUDA_MEM_INIT(d_seq1, seq1, strlen(seq1), char);
	
	char* d_mutant = NULL; // allocate mutant memory
	CUDA_MEM_INIT(d_mutant, mutant, mutant_size, char);
	
	float* d_weights = NULL; // allocate weights memory
	CUDA_MEM_INIT(d_weights, weights, WEIGHTS_NUM, float);
	
	float* d_similarity = NULL; // allocate similarity memory
	cudaMalloc((void**)&d_similarity, mutant_size * sizeof(float));
	
	int threads = (int)ceil(mutant_size);
	if (threads % 32 != 0) // assert threads is multiple of 32
		threads = threads+32 - threads%32; 
	
	int blocks = (mutant_size + threads - 1) / threads;
	
	calcSimilarityKernel<<<blocks, threads>>>(d_similarity, d_seq1, d_mutant, d_weights, mutant_size);

	cudaMemcpy(similarity, d_similarity, mutant_size * sizeof(float), cudaMemcpyDeviceToHost); // copy result to host memory
	
	cudaFree(d_seq1);
	cudaFree(d_mutant);
	cudaFree(d_weights);
	cudaFree(d_similarity);
}

// ========================

__device__ void getNK(int mutant_num, int seq2_len, int* n, int* k)
{
	int i;
	int num_of_mutants_in_row = seq2_len;

	for(i = 1; i < seq2_len; i++){
		if(mutant_num - (num_of_mutants_in_row - 1) > 0){
		    mutant_num -= (num_of_mutants_in_row - 1);
		    num_of_mutants_in_row--;
		}else{
		    break;
		}
	}
	
	*n = i;	
	*k = i + mutant_num;
}

__device__ float calcMutantScore(char* seq1, char* seq2, float* weights, float* d_conservative_matrix,int len2, int n, int k, int index, int offset)
{
	float score = 0;
	int i = 0, j = i;
	for (i = 0; i < len2 - 2; i++, j++)
	{
		if (j == n || j == k) 
			j++;
		float tmp_score = d_conservative_matrix[(seq1[i] - 'A') * 26 + (seq2[j] - 'A')];
		score += tmp_score;
	}	

	return score;	
}

__global__ void calcMutantBestScoreKernel(char* d_seq1, char* d_seq2, float* d_weights, float* d_mutantsBestScores, int* d_mutantsBestOffsets, int num_mutants, int maxOffset, int len2, float* d_conservative_matrix)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int n,k;
	float bestScore = -10000;
	int offset = 0;
	getNK(i+1, len2, &n, &k);
	if (i < num_mutants)
	{
		for (int j = 0; j < maxOffset; j++)
		{
			float score = calcMutantScore(&d_seq1[j], d_seq2, d_weights, d_conservative_matrix, len2, n, k, i, j);
			if (score > bestScore)
			{
				bestScore = score;
				offset = j;	
			}
		}
		d_mutantsBestScores[i]	= bestScore;
		d_mutantsBestOffsets[i] = offset;
	}
}

__global__ void reduction_cuda(float* d_mutantsBestScores, int* d_mutantsBestOffsets, float* d_reductionBestScores, int* d_reductionBestOffsets, int* d_reductionMutantNum, int num_of_elements){
	__shared__ int shared_mutant_num[256];
	__shared__ float shared_mutant_score[256];
	__shared__ float shared_mutant_best_offset[256];
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int tid = threadIdx.x;
	
	// load shared mem from global  mem
	shared_mutant_score[tid] = d_mutantsBestScores[i];
	shared_mutant_num[tid] = i;
	if(i >= num_of_elements)
		shared_mutant_score[tid] = 0;
	__syncthreads();	// make sure entine block is loaded!
	// do reduction in shared mem
	for(int s = blockDim.x / 2; s > 0; s >>= 1){
		if(tid < s){
			if(shared_mutant_score[tid] < shared_mutant_score[tid + s]){
				shared_mutant_score[tid] = shared_mutant_score[tid + s];
				shared_mutant_num[tid] = shared_mutant_num[tid + s];
				shared_mutant_best_offset[tid] = shared_mutant_best_offset[tid + s];
			}
		}
		__syncthreads();	// make sure all adds at ine stage are done!
	}
	// only thread 0 writes result for this block back to global mem
	if(tid == 0){
		d_reductionBestScores[blockIdx.x] = shared_mutant_score[tid];
		d_reductionBestOffsets[blockIdx.x] = shared_mutant_best_offset[tid];
		d_reductionMutantNum[blockIdx.x] = shared_mutant_num[tid];
	}
}

void calcBestScoreCUDA(char* seq1, char* seq2, float* conservative_matrix, float* weights)
{
	cudaError_t err = cudaSuccess;
	
	int len1 = strlen(seq1);
	int len2 = strlen(seq2);
	int maxOffset = len1 - (len2-2) + 1;
	int num_mutants = len2 * (len2 - 1) / 2;
	
	float* mutantsBestScores = (float*) malloc(num_mutants * sizeof(float));	
	int* mutantsBestOffsets = (int*) malloc(num_mutants * sizeof(int));		
	
	char* d_seq1 = NULL; // allocate seq1 memory
	CUDA_MEM_INIT(d_seq1, seq1, len1, char);
	
	char* d_seq2 = NULL; // allocate seq2 memory
	CUDA_MEM_INIT(d_seq2, seq2, len2, char);
	
	float* d_weights = NULL; // allocate weights memory
	CUDA_MEM_INIT(d_weights, weights, WEIGHTS_NUM, float);
	
	float* d_conservative_matrix = NULL; // allocate weights memory
	CUDA_MEM_INIT(d_conservative_matrix, conservative_matrix, 26 * 26, float);
	
	float* d_mutantsBestScores = NULL; 
	err = cudaMalloc((void**)&d_mutantsBestScores, num_mutants*sizeof(float));
	CUDA_ERR_CHECK(err, "Failed to allocate device memory");
	
	int* d_mutantsBestOffsets = NULL;
	err = cudaMalloc((void**)&d_mutantsBestOffsets, num_mutants*sizeof(int));
	CUDA_ERR_CHECK(err, "Failed to allocate device memory");
	
	int threads = 256;
	int blocks = (num_mutants + threads-1) / threads;
	calcMutantBestScoreKernel<<<blocks, threads>>>(d_seq1, d_seq2, d_weights, d_mutantsBestScores, d_mutantsBestOffsets, num_mutants, maxOffset, len2, d_conservative_matrix);
	
	cudaMemcpy(mutantsBestScores, d_mutantsBestScores, num_mutants * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(mutantsBestOffsets, d_mutantsBestOffsets, num_mutants * sizeof(int), cudaMemcpyDeviceToHost);
	
	/*
	float maxScore = -10000;
	int bestOffset = 0;
	int bestMutantNum = -1;
	
	for (int i = 0; i < num_mutants; i++)
	{
		if (mutantsBestScores[i] > maxScore)
		{
			maxScore = mutantsBestScores[i];
			bestOffset = mutantsBestOffsets[i];
			bestMutantNum = i;
		}
	}
	int n,k;
	getNK_CPU(bestMutantNum + 1, len2, &n, &k);
	printf("mutant num: %d, MS(%d,%d), score: %1.2f, offset: %d\n", bestMutantNum, n, k, maxScore, bestOffset);
	*/
	calcBestScoreOmp(mutantsBestScores, mutantsBestOffsets, num_mutants);
	
	free(mutantsBestScores);
	free(mutantsBestOffsets);
	cudaFree(d_seq1);
	cudaFree(d_seq2);
	cudaFree(d_mutantsBestScores);
	cudaFree(d_mutantsBestOffsets);
	cudaFree(d_conservative_matrix);
	
}
