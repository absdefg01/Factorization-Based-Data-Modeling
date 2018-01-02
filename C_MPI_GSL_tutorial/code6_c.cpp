#include <stdio.h>
#include <stdlib.h> 


struct myArray
{
  int len;
  double * data;
};


//在对应地址位置上的myArray上进行对应的操作 而不是建立一个新的
void fillArr(myArray * arr)
{
  for (int i = 0; i < arr->len; ++i)
  {
    arr->data[i] = i; //attention: "->" instead of "."
  }

}

//将建立一个新的myarray 进行操作
void fillArr2(myArray arr)
{
  for (int i = 0; i < arr.len; ++i)
  {
    arr.data[i] = i; 
  }

}

void printArr(myArray * arr){
  printf("%s\n", "myArray");
  for(int i = 0; i < arr->len; i++){
    printf("%f\n", arr->data[i]);
  }
}

int main(int argc, char** argv) 
{
   
  myArray arr;
  arr.len = 10;
  arr.data = (double *)malloc (arr.len * sizeof(double));

  fillArr(&arr);
  
  //how to print?
  printArr(&arr);
  
  return 0;
}







// printArr(&arr);

// void printArr(myArray * arr)
// {
//   for (int i = 0; i < arr->len; ++i)
//   {
//     printf("%f\n", arr->data[i]); 
//   }

// }