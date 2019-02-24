#ifndef __FUNC_EVAL_H_INCLUDED__
#define __FUNC_EVAL_H_INCLUDED__

#include <mpi.h>
#include <stdlib.h>
#include <math.h>

int InvertMatrix(int n, double *a, double *b, double *x1, double *x2, int taskid, int numtasks);

#endif /* __FUNC_EVAL_H_INCLUDED__ */
