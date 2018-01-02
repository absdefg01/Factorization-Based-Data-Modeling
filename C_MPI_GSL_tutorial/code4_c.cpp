#include <stdio.h>

//We need to include this header
#include <stdlib.h> 


int main(int argc, char** argv) {
   

  //This is an array 
  double myArray [10];

  for (int i = 0; i < 10; ++i)
  {
  	myArray[i] = i/2; //what happens here?
  }

  printf("myArray:\n");
  for (int i = 0; i < 10; ++i)
  {
  	printf("%f\n", myArray[i]);
  }

  //This is a dynamic array
  if(argc>1)
  {
  	int length = atoi(argv[1]); //character to integer conversion

  	double * myDynArray; // This is a pointer
  	
    //the length is int(argc)
  	myDynArray = (double *)malloc(length * sizeof(double)); //we allocate some memory

  	for (int i = 0; i < length; ++i)
  	{
  		myDynArray[i] = i;
  	}


  	printf("myDynArray:\n");
	for (int i = 0; i < length; ++i)
	{
		printf("%f\n", myDynArray[i]);
	}

  	free(myDynArray); //we must clean up the memory if we manually allocate


  }


  return 0;
}
