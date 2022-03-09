#pragma once

#include "mpi.h"
#include "header.h"

// ========================

void InitMPI(int* argc, char** argv[], int* rank, int* num_procs, MPI_Datatype* resultType, MPI_Datatype* bundleType);
void createResultType(MPI_Datatype* resultType);
void createBundleType(MPI_Datatype* bundleType);
void populateBundleData(char* seq1, char* seq2, float weights[], Bundle* bundle);

// ========================
