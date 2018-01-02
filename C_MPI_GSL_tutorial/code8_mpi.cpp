#include <mpi.h>
#include <stdio.h>
 #include <iostream>

int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
 
  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
 
 
  if(world_rank == 0)
  {
    //int length = ;
    //int * arr = (int*)malloc(length*sizeof(int));

    int a [] = {42,43,44};
    int destination = 1;

    //int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
    //buffer    initial address of send buffer (choice) 
    //count     number of elements in send buffer (nonnegative integer) 
    //datatype  datatype of each send buffer element (handle) 
    //dest      rank of destination (integer) 
    //tag       message tag (integer) 
    //comm      communicator (handle) 
    
    //&a : adress of send buffer
    //3 : length of buffer a
    //MPI_INT : datatype of each send buffer element
    //destination : rank of destination
    //0 : tag
    //MPI_COMM_WORLD : communicator
    MPI_Send(&a,3,MPI_INT,destination,0,MPI_COMM_WORLD);
    printf("Successfully sent!\n");
  }
  else
  {
    int b[3];
    MPI_Status st;
    int source = 0;
    //int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)
    //Output Parameters
    //buf    initial address of receive buffer (choice) 
    //status    status object (Status) 
    //Input Parameters
    //count    maximum number of elements in receive buffer (integer) 
    //datatype    datatype of each receive buffer element (handle) 
    //source    rank of source (integer) 
    //tag    message tag (integer) 
    //comm    communicator (handle) 

    //&b : adress of receive buffer
    //3 : length of receive buffer b
    //MPI_INT : datatype of each receive buffer element
    //source : rank of source
    //0 : tag
    //MPI_COMM_WORLD : communicator
    MPI_Recv(&b,3,MPI_INT,source,0,MPI_COMM_WORLD,&st);
    printf("Successfully received!\n");

    for(int i = 0; i<3; i++)
      printf("%d\n",b[i]); 


  }
 
  // Finalize the MPI environment.
  MPI_Finalize();
}
