#include <stdio.h>
 

int main(int argc, char** argv) {
   
  //The first argument is always the name of the executable
  // so argc = 1 + number of params entried by user.
  // %d : integer
  // %s : string
  // %f : float
  //etc ......
  printf("There are %d arguments.\n",argc-1);

  if(argc>1)
  {
  	printf("The arguments are:\n");
  	for (int i = 0; i < argc-1; ++i)
  	{
  		printf("\t%s\n",argv[i+1]);
  	}

  }

  return 0;
}
