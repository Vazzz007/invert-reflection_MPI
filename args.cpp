#include "args.h"

struct globalArgs_t {
    char *inFileName;           /* -i option */
    int matrixSize;             /* -n option */
    int verbose;                /* -v option */
    char *formula;              /* -f option */
    int max_out;                /* -m option */
} globalArgs;

static const char *optString = "i:n:f:vhm:?";

int get_args(int *matrixSize, char **inFileName, int *verbose, char **formula, 
             int *max_out, int argc, char **argv, int taskid
            ){
    
    int opt = 0;
    globalArgs.matrixSize = 0;
    globalArgs.inFileName = NULL;
    globalArgs.verbose = 0;
    globalArgs.formula = NULL;
    globalArgs.max_out = 5;
    opterr = 0;
    
    while( (opt = getopt( argc, argv, optString )) != -1 ) {
        switch( opt ) {
            case 'i':
                if ((globalArgs.matrixSize == 0) && (globalArgs.formula == NULL)) {
                    globalArgs.inFileName = optarg;
                } else {
                    if (taskid == 0) printf ("\nError: Only one of -i or -n/f can be used!\n");
                    if (taskid == 0) help();
                    return -1;
                }
                break;
            
            case 'n':
                if (globalArgs.inFileName == NULL) {
                    globalArgs.matrixSize = int(strtol(optarg, NULL, 10));
                } else {
                    if (taskid == 0) printf ("\nError: Only one of -i or -n/f can be used!\n");
                    if (taskid == 0) help();
                    return -1;
                }
                break;
            
            case 'v':
                globalArgs.verbose = 1;
                break;

            case 'f':
                if (globalArgs.inFileName == NULL) {
                    globalArgs.formula = optarg;
                } else {
                    if (taskid == 0) printf ("\nError: Only one of -i or -n/f can be used!\n");
                    if (taskid == 0) help();
                    return -1;
                }
                break;

            case 'h':
                if (taskid == 0) if (taskid == 0) help();
                return -1;
                break;

            case 'm':
                globalArgs.max_out = int(strtol(optarg, NULL, 10));
                break;

            case '?':
                if ((optopt == 'i') || (optopt == 'n') || (optopt == 'f') || (optopt == 'm')){
                    if (taskid == 0) printf ("\nError: Option -%c requires an argument.\n", optopt);
                }
                else if (isprint (optopt)){
                    if (taskid == 0) printf ("\nError: Unknown option `-%c'.\n", optopt);
                }
                else{
                    if (taskid == 0) printf (
                            "\nError: Unknown option character `\\x%x'.\n",
                            optopt);
                    if (taskid == 0) help();
                }      
                return -1;
            
            default:
                abort ();
                break;
        }
    }
    
    if ((globalArgs.inFileName == NULL) && ((globalArgs.matrixSize == 0) || 
	(globalArgs.formula == NULL))){
        if (taskid == 0) printf ("\nError: You must enter at least one of -i or -n/f\n");
        if (taskid == 0) help();
        return -1;
    }
    
    if ((globalArgs.matrixSize < 1) & (globalArgs.inFileName == NULL)){
        if (taskid == 0) printf ("\nError: Invalid matrix dimension!\n");
        return -1;
    }
    
    if (globalArgs.max_out < 1){
        if (taskid == 0) printf ("\nError: Invalid maximum output size!\n");
        return -1;
    }
    
    if (taskid == 0) printf ("\n\n inFileName = %s\n matrixSize = %d\n verbose = %d\n formula = %s\n",
              globalArgs.inFileName, globalArgs.matrixSize, globalArgs.verbose, globalArgs.formula);
    
    *matrixSize = globalArgs.matrixSize;
    *inFileName = globalArgs.inFileName;
    *verbose = globalArgs.verbose;
    *formula = globalArgs.formula;
    *max_out = globalArgs.max_out;
    
    return 0;
}

void help(){
    printf("\n\n/* Program supports the following command-line arguments:\n\
           * \n\
           * -i input_file_name.txt - name of the input file\n\
           * -n number - number of elements (default = 10)\n\
           * -v - option for debugging\n\
           * -f formula - define formula (choose from (choose from { sym ; smn ; glb ; 9 }\n\
           * -m number - maximum output size (default = 5)\n\n");
}
