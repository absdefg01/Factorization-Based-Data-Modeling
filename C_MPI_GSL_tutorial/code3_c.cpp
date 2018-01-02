#include <stdio.h>
 

//This is a function that doesn't return anything
void myFun1(int a)
{

  printf("The argument is %d\n",a);

}


//This is a function that returns something
double myFun2(int a)
{
  double b = a *1.5;
  return b;
}

int main(int argc, char** argv) {
   
  int a = 5;

  myFun1(a);
  double c = myFun2(a);
  printf("The output of fun 2 is %f\n", c);

  return 0;
}
