#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "create_matrix.h"
#include "func_eval.h"

static int max(int i, int j){
    if (i >= j) {
	return i;
    }
    return j;
}

static double sym_f(int i, int j){
    return (double)fabs(i - j);
}

static double symnul_f(int i, int j){
    return (double)fabs(i - j) + 1.0;
}

static double gilb_f(int i, int j){
    return 1.0 / (i + j + 1.0);
}

static double uniform(int i, int j, int n){
    return (double)(n - max(i, j));
}

int create_matrix(double *A, int n, char *formula, int taskid, int numtasks){
    
    int i, j;   
    int rows;

    if (taskid + 1 > n%numtasks) rows = n/numtasks;
    else rows = n/numtasks + 1; 
        
    if (strcmp(formula, "sym") == 0){
        for (i = 0; i < rows; ++i)
            for (j = 0; j < n; ++j)
                A[i * n + j] = sym_f(taskid + numtasks * i, j);
    } else if (strcmp(formula, "smn") == 0){
        for (i = 0; i < rows; ++i)
            for (j = 0; j < n; ++j)
                A[i * n + j] = symnul_f(taskid + numtasks * i, j);
    } else if (strcmp(formula, "glb") == 0){
        for (i = 0; i < rows; ++i)
            for (j = 0; j < n; ++j)
                A[i * n + j] = gilb_f(taskid + numtasks * i, j);
    } else if (strcmp(formula, "9") == 0){
        for (i = 0; i < rows; ++i)
            for (j = 0; j < n; ++j)
                A[i * n + j] = uniform(taskid + numtasks * i, j, n);
    } else {
        if (taskid == 0) printf ("Error: Invalid formula!\n");
        if (taskid == 0) help();
        //if (A != NULL) free(A);
        return -1;
    }
    
    return 0;
}

void PrintMatrix(int n, double *a, double *x, int max_out, int taskid, int numtasks)
{
    int i, j, m;
    MPI_Status status;

    m = (n < max_out) ? n : max_out;
    
    for (i = 0; i < m; i++)
    {

        if (taskid == 0)
        {
            if (taskid == i%numtasks)
            {
                printf("| ");
                for (j = 0; j < m; j++)
                    printf("%10.3e ", a[i/numtasks * n + j]);
                printf("|\n");
            }
            else
            {
                MPI_Recv(x, m, MPI_DOUBLE, i%numtasks, 0, MPI_COMM_WORLD, &status);
                printf("| ");
                for (j = 0; j < m; j++)
                    printf("%10.3e ", x[j]);
                printf("|\n");
            }
        }
        else if (taskid == i%numtasks)
        {
            for (j = 0; j < m; j++)
                x[j] = a[i/numtasks * n + j];
            MPI_Send(x, m, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
}

int InputMatrix(int n, double *A, FILE *fin, int taskid, int numtasks){
    int i, j, k;
    int rows;
    int sum = 0;
    double *aux = NULL;
    MPI_Status status;

    if (taskid == 0){

        if (taskid + 1 > n%numtasks) rows = n/numtasks;
        else rows = n/numtasks + 1;

        aux = (double *)malloc(n*n*sizeof(double));

        for (i = 0; i < n; ++i)
            for (j = 0; j < n; ++j)
                if (fscanf(fin, "%lf", &aux[i * n + j]) != 1)
                    return -1;

        for (i = 0; i < rows; ++i)
            for (j = 0; j < n; ++j)
                A[i * n + j] = aux[(taskid + numtasks * i)*n + j];
            //printf("\ntaskid = %d, rows = %d\n", taskid, rows);
            sum += rows;
            if (numtasks == 1) return 0;
            for (k = 1; k < numtasks; ++k){
                if (k + 1 > n%numtasks) rows = n/numtasks;
                else rows = n/numtasks + 1;
                //printf("\ntaskid %d waiting for sending %d elements(%f) to taskid %d\n", taskid, n*rows, aux[n*sum], k);
                for (int s = 0; s < rows; ++s){
                    MPI_Send(&aux[(k + numtasks * s)*n], n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD); // to WORKERS
                }
            
                //printf("\ntaskid %d sent %d elements(%f) to taskid %d\n", taskid, n*rows, aux[n*sum], k);
                sum += rows;
            }
    } else {
        
        if (taskid + 1 > n%numtasks) rows = n/numtasks;
        else rows = n/numtasks + 1;
        //printf("\ntaskid %d waiting to recieve %d elements from taskid 0", taskid, rows*n);
        for (int s = 0; s < rows; ++s){
            MPI_Recv(&A[s*n], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // from MASTER
        }
        //printf("\ntaskid %d recieved %d elements from taskid 0", taskid, rows*n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    free(aux);
    return 0;
}

double Residual(int n, double* a, double* b, int taskid, int numtasks)
{   
    int rows;
    double tmp;

    if (taskid + 1 > n%numtasks) rows = n/numtasks;
    else rows = n/numtasks + 1;

    double res = 0.0;

    for (int i = 0; i < rows; ++i){
        for (int j = 0; j < rows; ++j){

            tmp = 0.0;
            for (int k = 0; k < n; ++k)
                tmp += a[i * n + k] * b[j * n + k];
            //printf("tmp = %f\n", tmp);
            res += tmp * tmp;
            //printf("taskid = %d, i = %d\n", taskid, i);
        }
    //printf("\nmy_rank = %d, i = %d, residual = %10.3e\n", my_rank, i, *residual);
    }
    //printf("res = %f\n", res);
    
    return res;
}

void Transpose(int n, double* a, double* x, int taskid, int numtasks){
    int i, j, rows, sum = 0;
    double *aux = NULL;
    if (taskid == 0) aux = (double *)malloc(n*n*sizeof(double));
    MPI_Status status;
    
    for (i = 0; i < n; i++)
    {
        
        if (taskid == 0)
        {
            
            if (taskid == i%numtasks)
            {
                for (j = 0; j < n; j++)
                    aux[i*n + j] = a[i/numtasks * n + j];
            }
            else
            {
                MPI_Recv(x, n, MPI_DOUBLE, i%numtasks, 0, MPI_COMM_WORLD, &status);
                for (j = 0; j < n; j++){

                    aux[i*n + j] = x[j];
                    //printf(" %10.3e \n", aux[i*n + j]);
                }
            }
        }
        else if (taskid == i%numtasks)
        {
            for (j = 0; j < n; j++){
                //printf(" %10.3e \n", a[i/numtasks * n + j]);
                x[j] = a[i/numtasks * n + j];
            }
            MPI_Send(x, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }

        //printf("\ntaskid = %d, i = %d\n", taskid, i);
    }
    /*if(taskid == 0){
        for (i = 0; i < n; ++i){
            printf("| ");
            for (j = 0; j < n; ++j)
                printf("%10.3e ", aux[i * n + j]);
            printf("|\n");
        }
    }*/


    if (taskid + 1 > n%numtasks) rows = n/numtasks;
    else rows = n/numtasks + 1;
    if (taskid == 0){
        //printf("\ntaskid = %d\n", taskid);
        if (taskid + 1 > n%numtasks) rows = n/numtasks;
        else rows = n/numtasks + 1;

        for (i = 0; i < n; ++i)
            for (j = 0; j < rows; ++j)
                a[j * n + i] = aux[i*n + (j * numtasks + taskid)];
        //printf("\ntaskid = %d, rows = %d\n", taskid, rows);
        sum += rows;
        if (numtasks == 1) return;
        for (int k = 1; k < numtasks; ++k){
            if (k + 1 > n%numtasks) rows = n/numtasks;
            else rows = n/numtasks + 1;


            //printf("\ntaskid %d waiting for sending %d elements(%f) to taskid %d\n", taskid, n*rows, aux[n*sum], k);
            for (int s = 0; s < rows; ++s){
                for(i = 0; i < n; ++i){
                    x[i] = aux[i*n + (s * numtasks + k)];
                }
                MPI_Send(x, n, MPI_DOUBLE, k, 0, MPI_COMM_WORLD); // to WORKERS
            }
            
            //printf("\ntaskid %d sent %d elements(%f) to taskid %d\n", taskid, n*rows, aux[n*sum], k);
            sum += rows;
        }
    } else {
        
        if (taskid + 1 > n%numtasks) rows = n/numtasks;
        else rows = n/numtasks + 1;
        //printf("\ntaskid %d waiting to recieve %d elements from taskid 0", taskid, rows*n);
        for (int s = 0; s < rows; ++s){
            MPI_Recv(&a[s*n], n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status); // from MASTER
        }
        //printf("\ntaskid %d recieved %d elements from taskid 0", taskid, rows*n);
    }
    //printf("\ntaskid = %d\n", taskid);
    MPI_Barrier(MPI_COMM_WORLD);
    free(aux);
}