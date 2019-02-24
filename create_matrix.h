#ifndef __CREATE_MATRIX_H_INCLUDED__
#define __CREATE_MATRIX_H_INCLUDED__

#include <mpi.h>
#include <stdio.h>
#include <ctype.h>
#include <unistd.h>
#include <stdlib.h>
#include <fenv.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "args.h"

int create_matrix(double *A, int n, char *formula, int taskid, int numtasks);
void PrintMatrix(int n, double *a, double *x, int max_out, int taskid, int numtasks);
int InputMatrix(int n, double *A, FILE *fin, int taskid, int numtasks);
double SolutionError(int n, double* a, double* x);
int multi(int n, double* a, double* x, int my_rank, double *residual, int total_threads, double tmp);

#endif /* __CREATE_MATRIX_H_INCLUDED__ */
