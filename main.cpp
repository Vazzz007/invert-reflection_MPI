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
#include "args.h"
#include "create_matrix.h"
#include "func_eval.h"

#define EPS 1e-16

/* invert-reflection supports the following command-line arguments:
 * 
 * -i input_file_name.txt - name of the input file (default = NULL)
 * -n number - number of elements (default = 0)
 * -v - option for debugging
 * -f formula - define formula
 */

struct timespec diff(struct timespec start, struct timespec end);

struct timespec diff(struct timespec start, struct timespec end)
{
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

struct timespec time_thread_inv, time_thread_resid;

int main(int argc, char **argv){
    
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    
    int     n = 0,                          /* misc */
            rows,                           /* rows of matrix A sent to each worker */
            taskid,                         /* a task identifier */
            numtasks,                       /* number of tasks in partition */
            error_in = 0, error_out = 0,    /* error handlers */
            max_out = 0,                    /* size of output */
            verbose = 0,                    /* debug */
            temp = 0;
            

    double  *a,                             /* input matrix */
            *b,                             /* output matrix(inverse) */
            *x1, *x2;                       /* auxiliary vectors */
            //t,                              /* time */
            //res = 0.0;                      /* residual */
    
    char    *inFileName = NULL,             /* name of the input file */
            *formula = NULL;                /* unique name of matrix generator */

    FILE *fin;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    //MPI_Status status;

    formula = (char *)malloc(3 * sizeof(char));
    int inp_type = 0;
    
    if (taskid == 0)                        /* Master */
    {
        if ((get_args(&n, &inFileName, &verbose, &formula,
            &max_out, argc, argv, taskid, numtasks)) != 0
            )
        {
            if (taskid == 0) printf ("\nError: Can't get arguments!\n");
            MPI_Abort(MPI_COMM_WORLD, temp);
            return -1;
        }

        if (formula != NULL) inp_type = 0;
        else inp_type = 1;
    }

    MPI_Bcast(&inp_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (inp_type == 0) MPI_Bcast(formula, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (inp_type == 0) MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_out, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (taskid + 1 > n%numtasks) rows = n/numtasks;
    else rows = n/numtasks + 1;

    if (inp_type == 0){
        a = (double *)malloc (rows * n * sizeof(double));
        b = (double *)malloc (rows * n * sizeof(double));
        x1 = (double *)malloc (n * sizeof(double));
        x2 = (double *)malloc (n * sizeof(double));

        if ((a && b && x1 && x2) != 1) error_in = 1;

        MPI_Allreduce(&error_in, &error_out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        if (error_out != 0) {
            if (taskid == 0) printf ("\nError: Not enough memory!\n");
            if (a != NULL)
                free(a);
            if (b != NULL)
                free(b);
            if (x1 != NULL)
                free(x1);
            if (x2 != NULL)
                free(x2);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, temp);
            return -1;
        }
        
        if ((create_matrix(a, n, formula, taskid, numtasks)) != 0){
            if (taskid == 0) printf ("\nError: Can't create matrix!\n");
            if (a != NULL)
                free(a);
            if (b != NULL)
                free(b);
            if (x1 != NULL)
                free(x1);
            if (x2 != NULL)
                free(x2);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, temp);
            return -1;
        }
    } else {
        if (taskid == 0)                    /* Master */
        {
            fin = fopen(inFileName, "r");
            
            if (fin == NULL){
                if (taskid == 0) printf("\nError: Can't open file!\n");
                MPI_Abort(MPI_COMM_WORLD, temp);
                return -1;
            }
            
            if (fscanf(fin, "%d", &n) != 1){
                if (taskid == 0) printf("\nError: Can't read dimension from file!\n");
                fclose(fin);
                MPI_Abort(MPI_COMM_WORLD, temp);
                return -1;
            }
        
            if (n < 1){
                if (taskid == 0) printf ("\nError: Invalid matrix dimension!\n");
                fclose(fin);
                MPI_Abort(MPI_COMM_WORLD, temp);
                return -1;
            }
        }

        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (taskid + 1 > n%numtasks) rows = n/numtasks;
        else rows = n/numtasks + 1;

        a = (double *)malloc (rows * n * sizeof(double));
        b = (double *)malloc (rows * n * sizeof(double));
        x1 = (double *)malloc (n * sizeof(double));
        x2 = (double *)malloc (n * sizeof(double));

        if ((a && b && x1 && x2) != 1) error_in = 1;

        MPI_Allreduce(&error_in, &error_out, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
        if ((a && b && x1 && x2) != 1) {
            if (taskid == 0) printf ("\nError: Not enough memory!\n");
            if (a != NULL)
                free(a);
            if (b != NULL)
                free(b);
            if (x1 != NULL)
                free(x1);
            if (x2 != NULL)
                free(x2);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, temp);
            return -1;
        }

        /*if(taskid > 0){
            MPI_Send(a, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // to MASTER
        }*/

        if(taskid != -1){
            if (InputMatrix(n, a, fin, taskid, numtasks) != 0){
                if (taskid == 0) printf("\nError: Can't read matrix from file!\n");
                fclose(fin);
                free(a);
                free(b);
                free(x1);
                free(x2);
                MPI_Abort(MPI_COMM_WORLD, temp);
                return -1;
            }
        }

        /*if(taskid > 0) {
            if (taskid + 1 > n%numtasks) rows = n/numtasks;
            else rows = n/numtasks + 1;
            MPI_Recv(&a, rows*n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
        }*/
    }
    //printf("\ntaskid = %d\n", taskid);
    if (taskid == 0) printf("\nMatrix A:\n");
    
    PrintMatrix(n, a, x1, max_out, taskid, numtasks);

    if (taskid == 0) printf("\n");
    
    if (taskid == 0) printf("\nCalculating...\n");
    
/*    if( clock_gettime( CLOCK_MONOTONIC, &time_start) == -1 ) {
        perror( "clock gettime" );
        exit( EXIT_FAILURE );
    }
    */
// __________________Go into the main algorithm___________________
    MPI_Barrier(MPI_COMM_WORLD);

    InvertMatrix(n, a, b, x1, x2, taskid, numtasks);

    MPI_Barrier(MPI_COMM_WORLD);
// _____________________________Exit______________________________
        
    if (taskid == 0) printf("\nMatrix A^{-1}:\n");
    PrintMatrix(n, b, x1, max_out, taskid, numtasks);

    free(a);
    free(b);
    free(x1);
    free(x2);

    MPI_Finalize();

    return 0;
}


/*    if( clock_gettime( CLOCK_MONOTONIC, &time_end) == -1 ) {
        perror( "clock gettime" );
        exit( EXIT_FAILURE );
    }
    time_end = diff(time_start, time_end);
    time_total = time_end;
    
    if (result == -1){
        printf("\nCan't invert.\n");
    } else {
        printf("\nMatrix A^{-1}:\n");
        PrintMatrix(n, X, max_out);
        printf("\n");
        
        printf("\nInversion time \t\t= %f sec.\nInversion_thread_time\t= %f sec.\n\n",
               (double)time_total.tv_sec + (double)time_total.tv_nsec/(double)1000000000,
               (double)time_thread_inv.tv_sec + (double)time_thread_inv.tv_nsec/(double)1000000000);
        
        
        if (inFileName == NULL){
            if ((create_matrix(A, n, formula)) != 0){
                printf ("\nError: Can't create matrix!\n");
                free(A);
                free(X);
                free(threads);
                free(args);
                return -1;
            }
        } else {
            fin = fopen(inFileName, "r");

            if (fin == NULL){
                printf("\nError: Can't open file!\n");
                free(A);
                free(X);
                free(threads);
                free(args);
                return -1;
            }

            if (fscanf(fin, "%d", &n) != 1){
                printf("\nError: Can't read dimension from file!\n");
                free(A);
                free(X);
                free(threads);
                free(args);
                fclose(fin);
                return -1;
            }

            if (n < 1){
                printf ("\nError: Invalid matrix dimension!\n");
                free(A);
                free(X);
                free(threads);
                free(args);
                fclose(fin);
                return -1;
            }
        
            if (InputMatrix(n, A, fin) != 0){
                printf("\nError: Can't read matrix from file!\n");
                fclose(fin);
                free(A);
                free(X);
                free(threads);
                free(args);
                return -1;
            }
            
            fclose(fin);
        }

        args_mul = (ARGS_mul*)malloc(total_threads * sizeof(ARGS_mul));

        if( clock_gettime( CLOCK_MONOTONIC, &time_start) == -1 ) {
            perror( "clock gettime" );
            exit( EXIT_FAILURE );
        }
        
        for (int i = 0; i < total_threads; i++){
            args_mul[i].n = n;
            args_mul[i].A = A;
            args_mul[i].X = X;
            args_mul[i].my_rank = i;
            args_mul[i].residual = 0.0;
            args_mul[i].total_threads = total_threads;
            args_mul[i].tmp = 0.0;
        }

        for (int i = 0; i < total_threads; i++){
            if (pthread_create(threads + i, NULL, Multiplication, args_mul + i)){
                printf("Cannot create thread %d!\n", i);

                if (A) free(A);
                if (X) free(X);
                if (threads) free(threads);
                if (args) free(args);
                if (args_mul) free(args_mul);

                return -1;
            }
        }

        for (int i = 0; i < total_threads; i++){
            if (pthread_join(threads[i], NULL)){
                printf("Cannot wait thread %d!\n", i);

                if (A) free(A);
                if (X) free(X);
                if (threads) free(threads);
                if (args) free(args);
                if (args_mul) free(args_mul);

                return -1;
            }
        }

        if( clock_gettime( CLOCK_MONOTONIC, &time_end) == -1 ) {
            perror( "clock gettime" );
            exit( EXIT_FAILURE );
        }
        time_end = diff(time_start, time_end);
        printf("\nResidual time \t\t= %f sec.\nResidual_thread_time\t= %f sec.\n\n",
               (double)time_end.tv_sec + (double)time_end.tv_nsec/(double)1000000000,
               (double)time_thread_resid.tv_sec + (double)time_thread_resid.tv_nsec/(double)1000000000);
        time_total.tv_sec += time_end.tv_sec;
        time_total.tv_nsec += time_end.tv_nsec;

        time_thread_total.tv_sec = time_thread_inv.tv_sec + time_thread_resid.tv_sec;
        time_thread_total.tv_nsec = time_thread_inv.tv_nsec + time_thread_resid.tv_nsec;


        for (int i = 0; i < total_threads; i++){
            nev += (1) * args_mul[i].residual;
        }

        printf("\nTotal time \t\t= %f sec.\nTotal_thread_time\t= %f sec.\n\n",
               (double)time_total.tv_sec + (double)time_total.tv_nsec/(double)1000000000,
               (double)time_thread_total.tv_sec + (double)time_thread_total.tv_nsec/(double)1000000000);

        printf("\nSolution threaded ||A * A^{-1} - I||\t= %e\n", nev);

        free(args_mul);

        //printf("\nSolution ||A * A^{-1} - I||\t= %e\n", SolutionError(n, A, X));
    }

    free(threads);
    free(args);
    free(A);
    free(X);
    
    return 0;
}
*/