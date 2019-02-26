#include <math.h>
#include <mpi.h>

#include "func_eval.h"

int InvertMatrix(int n, double *a, double *b, double *x1, double *x2, int taskid, int numtasks)
{
  int i, j, k;
  int rows;
  int first;
  double tmp;
  double q1, q2;
  double inv_norma;

  if (taskid + 1 > n%numtasks) rows = n/numtasks;
  else rows = n/numtasks + 1;

  for (i = 0; i < rows; i++)
    for (j = 0; j < n; j++)
      b[i * n + j] = (taskid + i * numtasks == j);

  for (i = 0; i < n; i++)
  {
    if (taskid == i%numtasks)
    {
      if (i == n - 1)
      {
        //printf("\na = %f\n", a[i/numtasks * n + i]);
        if(fabs(a[i/numtasks * n + i]) < 1e-50){
          return -1;
        }

        tmp = 1.0/a[i/numtasks * n + i];
        a[i/numtasks * n + i] = 1.0;
        for (j = 0; j < n; j++)
          b[i/numtasks * n + j] *= tmp;
        continue;
      }

      for (q1 = 0.0, j = i/numtasks + 1; j < rows; j++)
        q1 += a[j * n + i] * a[j * n + i];

      MPI_Reduce(&q1, &q2, 1, MPI_DOUBLE, MPI_SUM, i%numtasks, MPI_COMM_WORLD);

      tmp = sqrt(q2 + a[i/numtasks * n + i] * a[i/numtasks * n + i]);
      a[i/numtasks * n + i] -= tmp;

      if(sqrt(q2 + a[i/numtasks * n + i] * a[i/numtasks * n + i]) < 1e-50){
          return -1;
        }
      inv_norma = 1.0/sqrt(q2 + a[i/numtasks * n + i] * a[i/numtasks * n + i]);
      MPI_Bcast(&inv_norma, 1, MPI_DOUBLE, i%numtasks, MPI_COMM_WORLD);

      for (j = i/numtasks; j < rows; j++)
        a[j * n + i] *= inv_norma;

      for (k = i + 1; k < n; k++)
      {
        for (q1 = 0.0, j = i/numtasks; j < rows; j++)
          q1 += a[j * n + i] * a[j * n + k];
        x1[k] = q1;
      }

      MPI_Allreduce(x1, x2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (k = i + 1; k < n; k++)
      {
        x2[k] *= 2.0;
        for (j = i/numtasks; j < rows; j++)
          a[j * n + k] -= x2[k] * a[j * n + i];
      }

      for (k = 0; k < n; k++)
      {
        for (q1 = 0.0, j = i/numtasks; j < rows; j++)
          q1 += a[j * n + i] * b[j * n + k];
        x1[k] = q1;
      }

      MPI_Allreduce(x1, x2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (k = 0; k < n; k++)
      {
        x2[k] *= 2.0;
        for (j = i/numtasks; j < rows; j++)
          b[j * n + k] -= x2[k] * a[j * n + i];
      }

      a[i/numtasks * n + i] = tmp;
      for (j = i/numtasks + 1; j < rows; j++)
        a[j * n + i] = 0.0;

      tmp = 1.0/a[i/numtasks * n + i];
      for (j = i; j < n; j++)
        a[i/numtasks * n + j] *= tmp;
      for (j = 0; j < n; j++)
        b[i/numtasks * n + j] *= tmp;
    }
    else
    {
      if (i == n - 1) continue;

      if (taskid > i%numtasks) first = i/numtasks;
      else first = i/numtasks + 1;

      for (q1 = 0.0, j = first; j < rows; j++)
        q1 += a[j * n + i] * a[j * n + i];

      MPI_Reduce(&q1, &q2, 1, MPI_DOUBLE, MPI_SUM, i%numtasks, MPI_COMM_WORLD);

      MPI_Bcast(&inv_norma, 1, MPI_DOUBLE, i%numtasks, MPI_COMM_WORLD);

      for (j = first; j < rows; j++)
        a[j * n + i] *= inv_norma;

      for (k = i + 1; k < n; k++)
      {
        for (q1 = 0.0, j = first; j < rows; j++)
          q1 += a[j * n + i] * a[j * n + k];
        x1[k] = q1;
      }

      MPI_Allreduce(x1, x2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (k = i + 1; k < n; k++)
      {
        x2[k] *= 2.0;
        for (j = first; j < rows; j++)
          a[j * n + k] -= x2[k] * a[j * n + i];
      }

      for (k = 0; k < n; k++)
      {
        for (q1 = 0.0, j = first; j < rows; j++)
          q1 += a[j * n + i] * b[j * n + k];
        x1[k] = q1;
      }

      MPI_Allreduce(x1, x2, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      for (k = 0; k < n; k++)
      {
        x2[k] *= 2.0;
        for (j = first; j < rows; j++)
          b[j * n + k] -= x2[k] * a[j * n + i];
      }

      for (j = first; j < rows; j++)
        a[j * n + i] = 0.0;
    }
  }

  for (i = n - 1; i >= 1; i--)
  {
    if (taskid == i%numtasks)
    {
      MPI_Bcast(b + i/numtasks * n, n, MPI_DOUBLE, i%numtasks, MPI_COMM_WORLD);
      for (j = i/numtasks - 1; j >= 0 ; j--)
        for (k = 0; k < n; k++)
          b[j * n + k] -= b[i/numtasks * n + k] * a[j * n + i];
    }
    else
    {
      MPI_Bcast(x1, n, MPI_DOUBLE, i%numtasks, MPI_COMM_WORLD);

      if (taskid < i%numtasks) first = i/numtasks;
      else if (i/numtasks - 1 >= 0) first = i/numtasks - 1;
      else continue;

      for (j = first; j >= 0; j--)
        for (k = 0; k < n; k++)
          b[j * n + k] -= x1[k] * a[j * n + i];
    }
  }

  return 0;
}
