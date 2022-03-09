#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "omp.h"
// ===================== macros

#define NUM_PROCS 2

#define MASTER 0

#define READ_FILE_NAME "input.txt"
#define WRITE_FILE_NAME "output.txt"

#define SEQ1_LEN 5000
#define SEQ2_LEN 3000

#define WEIGHTS_NUM 4
#define ABC_NUM 26

#define CONSERVATIVES_LEN 9
#define SEMI_CONSERVATIVES_LEN 11

// ===================== 

typedef struct
{
	float score;
	int offset, n, k;

} Result;

typedef struct
{
	char seq1[SEQ1_LEN];
	int seq1_len;
	char seq2[SEQ2_LEN];
	int seq2_len;
	float weights[WEIGHTS_NUM];
	
} Bundle;

// ===================== fileUtil

char** readFromFile(const char* fileName, float weights[], char* seq1, int* numOfSeq2);
void writeToFile(const char* fileName, Result* results, int num_seqs);

// ===================== cFunctions

Result calcBestScore(Bundle bundle, int rank);
Result calcOffsetBestResult(char* seq1, char* seq2, float weights[], int num_mutants, int offset);
Result calcMutantResultOMP(float* similarity, int mutant_size, int offset, int n, int k);

// ===================== helpers

void assignOffsets(int* start, int* end, int maxOffset, int rank);
char* MS(char* seq, int n, int k);
int** createNKs(int num_mutants);

// ===================== cudaFunctions

void calcMutantSimilarityCUDA(float* similarity, char* seq1, char* mutant, float weights[]);
void calcBestScoreCUDA(char* seq1, char* seq2, float* conservative_matrix, float* weights);
void printBundle(Bundle bundle);


void calcBestScoreOmp(float* mutantsBestScores, int* mutantsBestOffsets, int num_mutants);

