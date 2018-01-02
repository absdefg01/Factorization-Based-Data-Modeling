#include <stdio.h>
#include <stdlib.h> 


void printArray (double * arr, int length)
{

  for (int i = 0; i < length; ++i)
  {
    printf("%f\n", arr[i]);
  }

}

//printArray(&(myArray[0]),10);


int main(int argc, char** argv) {
   

  //This is an array 
  double myArray [10];

  for (int i = 0; i < 10; ++i)
  {
  	myArray[i] = i/2; //what happens here?
  }

  printf("myArray:\n");
  printArray (&(myArray[0]), 10);

  //This is a dynamic array
  if(argc>1)
  {
  	int length = atoi(argv[1]); //character to integer conversion

  	double * myDynArray; // This is a pointer
  	
  	myDynArray = (double *)malloc(length * sizeof(double)); //we allocate some memory

  	for (int i = 0; i < length; ++i)
  	{
  		myDynArray[i] = i;
  	}


  	printf("myDynArray:\n");
	  printArray (myDynArray, length);

  	free(myDynArray); //we must clean up the memory if we manually allocate


  }


  return 0;
}
