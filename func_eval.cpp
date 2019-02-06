#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "mpi.h"

#include "func_eval.h"

#define MASTER 0               /* taskid of first task */
#define FROM_MASTER 1          /* setting a message type */
#define FROM_WORKER 2          /* setting a message type */

int InvMatrix(int n, double *a, double *x, int total_threads, int *result, int v)
{
// ________Declaration and initialization of mpi parameters_____________

int	numtasks,              /* number of tasks in partition */
	taskid,                /* a task identifier */
	numworkers,            /* number of worker tasks */
	source,                /* task id of message source */
	dest,                  /* task id of message destination */
	mtype,                 /* message type */
	rows,                  /* rows of matrix A sent to each worker */
	averow, extra, offset, /* used to determine rows sent to each worker */
	i, j, k, rc;           /* misc */
	first_row;
	last_row;

double tmp1, tmp2;

MPI_Status status;

int argc = 1;
char argv[1][9] = {
	{'-', 'n', ' ', total_threads}
}

MPI_Init(&argc,&argv);
MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
if (numtasks < 2 ) {
	printf("Need at least two MPI tasks. Quitting...\n");
	MPI_Abort(MPI_COMM_WORLD, rc);
	exit(1);
	}
numworkers = numtasks-1;

// ________________End_______________________________


/**************************** master task open ************************************/

if (taskid == MASTER)
   {
      printf("mpi_invertion has started with %d tasks.\n",numtasks);

      /* Measure start time */
      double start = MPI_Wtime();

      for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				x[i * n + j] = (double)(i == j);

      /* Send matrix data to the worker tasks */
      averow = NRA/numworkers;
      extra = NRA%numworkers;
      offset = 0;
      mtype = FROM_MASTER;
      for (dest=1; dest<=numworkers; dest++)
      {
         rows = (dest <= extra) ? averow+1 : averow;   	
         if (v == 1)
         	printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
         MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
         MPI_Send(&a[offset][0], rows*n, MPI_DOUBLE, dest, mtype,
                   MPI_COMM_WORLD);
         MPI_Send(&b, n*n, MPI_DOUBLE, dest, mtype, MPI_COMM_WORLD);
         offset = offset + rows;
      }

      /* Receive results from worker tasks */
      mtype = FROM_WORKER;
      for (i=1; i<=numworkers; i++)
      {
         source = i;
         MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
         MPI_Recv(&c[offset][0], rows*n, MPI_DOUBLE, source, mtype, 
                  MPI_COMM_WORLD, &status);
         // printf("Received results from task %d\n",source);
      }

      /* Print results */
      /*
      printf("******************************************************\n");
      printf("Result Matrix:\n");
      for (i=0; i<NRA; i++)
      {
         printf("\n"); 
         for (j=0; j<n; j++) 
            printf("%6.2f   ", c[i][j]);
      }
      printf("\n******************************************************\n");
      */

      /* Measure finish time */
      double finish = MPI_Wtime();
      printf("Done in %f seconds.\n", finish - start);
   }

/**************************** worker task ************************************/
if (taskid > MASTER)
   {
      mtype = FROM_MASTER;
      MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&a, rows*NCA, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&b, NCA*NCB, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);

      for (k=0; k<NCB; k++)
         for (i=0; i<rows; i++)
         {
            c[i][k] = 0.0;
            for (j=0; j<NCA; j++)
               c[i][k] = c[i][k] + a[i][j] * b[j][k];
         }
      mtype = FROM_WORKER;
      MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
      MPI_Send(&c, rows*NCB, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
   }



	for (i = 0; i < n-1; i++){

        if (v == 1)
            printf("\nmy_rank = %d , i = %d \n", taskid, i);
        //Нужно научиться передавать данные в другой thread
        synchronize(total_threads);
		tmp1 = 0.0;
		for (j = i + 1; j < n; j++)
			tmp1 += a[j * n + i] * a[j * n + i];

		tmp2 = sqrt(tmp1 + a[i * n + i] * a[i * n + i]);
            
            
        if (tmp2 < 1e-100){
            *result = -1;
            //printf ("status in func = %d", *status);
            pthread_exit(NULL);
            return -1;
        }
        synchronize(total_threads);
        if (taskid == 1){
			a[i * n + i] -= tmp2;
        }
        synchronize(total_threads);
        tmp1 = sqrt(tmp1 + a[i * n + i] * a[i * n + i]);

        synchronize(total_threads);
        if (v == 1)
            printf("\nmy_rank = %d , tmp1 = %f \n", my_rank, tmp1);

        if (tmp1 < 1e-100){
            if (taskid == 1)
                a[i * n + i] += tmp2;
            //Здесь нужно передать индекс i другим процессам
            printf("Change index");
            continue;
        }

        synchronize(total_threads);

        tmp1 = 1.0/tmp1;
        if (my_rank == 0){
            
            for (j = i; j < n; j++)
                a[j * n + i] *= tmp1;
            
        }
        

		synchronize(total_threads);

		first_row = (n - i - 1) * my_rank;
		first_row = first_row/total_threads + i + 1;
		last_row = (n - i - 1) * (my_rank + 1);
		last_row = last_row/total_threads + i + 1;

		for (k = first_row; k < last_row; k++)
		{
			tmp1 = 0.0;
			for (j = i; j < n; j++)
				tmp1 += a[j * n + i] * a[j * n + k];

			tmp1 *= 2.0;
			for (j = i; j < n; j++)
				a[j * n + k] -= tmp1 * a[j * n + i];
		}
		synchronize(total_threads);

		first_row = n * my_rank;
		first_row = first_row/total_threads;
		last_row = n * (my_rank + 1);
		last_row = last_row/total_threads;

		for (k = first_row; k < last_row; k++)
		{
			tmp1 = 0.0;
			for (j = i; j < n; j++)
				tmp1 += a[j * n + i] * x[j * n + k];

			tmp1 *= 2.0;
			for (j = i; j < n; j++)
				x[j * n + k] -= tmp1 * a[j * n + i];

		}
		synchronize(total_threads);

		if (my_rank == 0)
			a[i * n + i] = tmp2;

	}
    synchronize(total_threads);

	first_row = n * my_rank;
	first_row = first_row/total_threads;
	last_row = n * (my_rank + 1);
	last_row = last_row/total_threads;
    if (v == 1)
        printf("\nmy_rank = %d , Go to gauss\n", my_rank);

	for (k = first_row; k < last_row; k++){
		for (i = n - 1; i >= 0; i--)
		{
			tmp1 = x[i * n + k];
			for (j = i + 1; j < n; j++)
				tmp1 -= a[i * n + j] * x[j * n + k];
			x[i * n + k] = tmp1/a[i * n + i];
		}
        if (v == 1)
            printf("\nmy_rank = %d , k = %d \n", my_rank, k);
    }

    if (v == 1)
        printf("\nmy_rank = %d , Work done!\n", my_rank);
	return 0;
	MPI_Finalize();
}