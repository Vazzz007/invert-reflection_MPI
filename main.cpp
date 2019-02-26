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
            *x1, *x2,                       /* auxiliary vectors */
            total, start, end,              /* time */
            res = 0.0, sumres;              /* residual */
    
    char    *inFileName = NULL,             /* name of the input file */
            *formula = NULL;                /* unique name of matrix generator */

    FILE *fin;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    //MPI_Status status;

    if (taskid > 0) formula = (char *)malloc(3 * sizeof(char));
    printf("\n1taskid = %d, formula = %p\n", taskid, formula);
    int inp_type = 0;
    
    if (taskid == 0)                        /* Master */
    {
        if ((get_args(&n, &inFileName, &verbose, &formula,
            &max_out, argc, argv, taskid)) != 0
            )
        {
            if (taskid == 0) printf ("\nError: Can't get arguments!\n");
            MPI_Abort(MPI_COMM_WORLD, temp);
            return -1;
        }

        if (formula != NULL) inp_type = 0;
        else inp_type = 1;
    }
    printf("\n2taskid = %d, formula = %p\n", taskid, formula);

    MPI_Bcast(&inp_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (inp_type == 0) MPI_Bcast(formula, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (inp_type == 0) MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_out, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("\n3taskid = %d, formula = %p\n", taskid, formula);

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

    //printf("\ntaskid = %d\n", taskid);
    if (taskid == 0) printf("\nMatrix A:\n");
    
    PrintMatrix(n, a, x1, max_out, taskid, numtasks);

    if (taskid == 0) printf("\n");
    
    if (taskid == 0) printf("\nCalculating...\n");

    if (numtasks > n) {
        if (taskid == 0) printf("\nSet less(than n = %d) processes!\n", n + 1);
        free(a);
        free(b);
        free(x1);
        free(x2);
        MPI_Abort(MPI_COMM_WORLD, temp);
        return -1;
    }


// __________________Go into the main algorithm___________________
    MPI_Barrier(MPI_COMM_WORLD);

    start = MPI_Wtime();

    error_in = InvertMatrix(n, a, b, x1, x2, taskid, numtasks);

    end = MPI_Wtime() - start;

    if(error_in != 0){
        printf("\nCan't invert.\n");
        free(a);
        free(b);
        free(x1);
        free(x2);
        MPI_Abort(MPI_COMM_WORLD, temp);
        return -1;
    }

    MPI_Barrier(MPI_COMM_WORLD);
// _____________________________Exit______________________________
        
    MPI_Allreduce(&end, &total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    if (0){
        if (taskid == 0) printf("\nCan't invert.\n");
    } else {
        if (taskid == 0) printf("\nMatrix A^{-1}:\n");
        PrintMatrix(n, b, x1, max_out, taskid, numtasks);
        if (taskid == 0) printf("\n");
        
        if (taskid == 0) printf("\nInversion time \t\t= %f sec.\nInversion_thread_time\t= %f sec.\n\n",
               (double)total,
               (double)end);
        
        
        if (inp_type == 0){
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

            if(taskid == 0) fclose(fin);
        }

        Transpose(n, b, x1, taskid, numtasks);

        //PrintMatrix(n, b, x1, max_out, taskid, numtasks);

        MPI_Barrier(MPI_COMM_WORLD);

        start = MPI_Wtime();

        res = Residual(n, a, b, taskid, numtasks);

        end = MPI_Wtime() - start;

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Allreduce(&res, &sumres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        MPI_Allreduce(&end, &total, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        //printf("res = %f\n", res);
        //printf("sumres = %f\n", sumres);

        if(taskid == 0) printf("\nResidual time \t\t= %f sec.\nResidual_thread_time\t= %f sec.\n\n",
               (double)total,
               (double)end);

        /*printf("\nTotal time \t\t= %f sec.\nTotal_thread_time\t= %f sec.\n\n",
               (double)time_total.tv_sec + (double)time_total.tv_nsec/(double)1000000000,
               (double)time_thread_total.tv_sec + (double)time_thread_total.tv_nsec/(double)1000000000);
*/
        if(taskid == 0) printf("\nSolution threaded ||A * A^{-1} - I||\t= %e\n", sqrt(fabs(sumres - n)));
    }

    free(x1);
    free(x2);
    free(a);
    free(b);

    printf("\n4taskid = %d, formula = %p\n", taskid, formula);
    //if (taskid > 0) free(formula);

    MPI_Finalize();
    
    return 0;
}
