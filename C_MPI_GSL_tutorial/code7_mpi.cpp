#include <mpi.h>
#include <stdio.h>
 

//to compile: mpicxx code7_mpi.cpp -o code7_mpi
//to run: mpirun -n "number of processes" ./code7_mpi

//world rank is numero of processor
//workd size is the number of processes
int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int world_size;

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  //here : we can notice that world_size = number of processes
  /*
  printf("%s\n", "world_size");
  printf("%d\n", world_size);
  */

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
 
  // Print off a hello world message
  printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, world_rank, world_size);
 
  /*
  if(world_rank == 0){
    printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, world_rank, world_size);
  }
  */

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}


/*
#include <mpi.h>
#include <stdio.h>
 

//to compile: mpicxx code7_mpi.cpp -o code7_mpi
//to run: mpirun -n "number of processes" ./code7_mpi

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  if(world_rank == 0){
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         world_rank, world_size);
  }else{
    // Print off a hello world message
    printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         world_rank, world_size);
  } 

  
  // Print off a hello world message
  printf("Hello world from processor %s, rank %d"
         " out of %d processors\n",
         processor_name, world_rank, world_size);
 

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}

*/
