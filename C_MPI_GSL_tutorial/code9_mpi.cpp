#include <mpi.h>
#include <stdio.h>
#include <stdlib.h> 
 #include <fstream>
 #include <iostream>
int main(int argc, char** argv) {
  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
 
  // Get the number of processes
  int numProcesses;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
 
  // Get the rank of the process
  int processId;
  MPI_Comm_rank(MPI_COMM_WORLD, &processId);
 
 
  int arr[3];
  int size = 3;
  for (int i = 0; i < 3; ++i)
  {
    arr[i] = processId*10 + i;
  }

  int * tempBuf = (int *)malloc(size*sizeof(int));

  int dest = (processId-1+numProcesses)%numProcesses;
  int src  = (processId+1)%numProcesses;


  printf("Process %d is sending %d, %d, %d to Process %d\n",processId, arr[0],arr[1],arr[2], dest);

  MPI_Status st;
  //int MPI_Sendrecv
  //const void *sendbuf, int sendcount, MPI_Datatype sendtype,
  //int dest, int sendtag,void *recvbuf, int recvcount, 
  //MPI_Datatype recvtype,int source, int recvtag,MPI_Comm comm, MPI_Status *status

  //sendbuf    initial address of send buffer (choice) 
  //sendcount    number of elements in send buffer (integer) 
  //sendtype    type of elements in send buffer (handle) 
  //dest    rank of destination (integer) 
  //sendtag    send tag (integer) 
  //recvcount    number of elements in receive buffer (integer) 
  //recvtype    type of elements in receive buffer (handle) 
  //source    rank of source (integer) 
  //recvtag    receive tag (integer) 
  //comm    communicator (handle) 

  MPI_Sendrecv(arr, size, MPI_INT,dest, 0,
              tempBuf, size, MPI_INT,src, 0,
              MPI_COMM_WORLD, &st);

  //copy incoming data
  memcpy(arr, tempBuf, size*sizeof(int));

  printf("Process %d is received %d, %d, %d from Process %d\n",processId, arr[0],arr[1],arr[2], src);

  // Finalize the MPI environment.
  MPI_Finalize();

  free(tempBuf);
}
