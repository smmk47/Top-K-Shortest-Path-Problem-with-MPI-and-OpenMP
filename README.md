 # MPI and OpenMP Code README
 ## Overview:

 This README provides instructions for compiling and executing an MPI (Message Passing Interface) and OpenMP (Open Multi 
  Processing) hybrid code on a Linux environment, for the PDC PROJECT of students
  
Saad Ahmed Qureshi 21i-0616

Murtaza Kazmi 21i-0685

Ibrahim Salman 21i-2516

## Prerequisites:
MPI and OpenMP libraries installed on your system. You can typically install them using your package manager. For example, on Ubuntu, you can install MPI and OpenMP libraries with:


#### sudo apt-get install mpich libomp-dev
## Compilation:
Navigate to the directory containing the MPI and OpenMP code files.
Open a terminal window.
Use the MPI compiler wrapper (mpicc) along with OpenMP flags to compile the code. The general syntax is:

#### mpicc -o <executable_name> <source_files> -fopenmp
### Example:
#### mpicc -o my_mpi_openmp_program mpi_openmp_code.c -fopenmp
Upon successful compilation, an executable file will be generated (in this example, my_mpi_openmp_program).
## Execution:
After compiling the code, you can execute the MPI and OpenMP program using the mpirun command.
Open a terminal window.
Use the following syntax to run the MPI and OpenMP program:

#### mpirun -np <num_processes> <executable_name>
Replace <num_processes> with the desired number of MPI processes you want to run. 
Replace <executable_name> with the name of the executable generated during compilation.
### Example:

#### mpirun -np 4 my_mpi_openmp_program.
This command will execute the MPI and OpenMP program using 4 MPI processes.
## Notes:
Ensure that the MPI and OpenMP program is designed to run in parallel and has appropriate MPI and OpenMP function calls for communication between processes and threads. 

Make sure to adjust the number of MPI processes (-np) based on your system's capabilities and the requirements of your MPI and OpenMP program.
