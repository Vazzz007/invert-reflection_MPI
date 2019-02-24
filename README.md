# invert-reflection_MPI
The MPI reflection method for finding the inverse matrix using an adjugate matrix.

Use "make" to compile.

Program supports the following command-line arguments:
  * -i input.txt - name of the input file(exactly)
  * -n number - number of elements (default = 10)
  * -v - option for debugging
  * -f formula - define formula (choose from { sym ; smn ; glb ; 9 }
  * -m number - maximum output size (default = 5)
  
  Debug: mpirun -np 2 xterm -e gdb --args ./MPI -n 200 -f sym -m 7
  
# The samples of using:
  
  make
  
  mpirun -np 4 ./MPI -n 200 -f sym -m 7
  
  mpirun -np 2 ./MPI -i input.txt
  
  make clean
