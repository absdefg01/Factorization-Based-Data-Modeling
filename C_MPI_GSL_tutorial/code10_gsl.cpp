#include <stdio.h>
#include <gsl/gsl_matrix.h>

//to compile: g++ code10_gsl.cpp -Wall -lgsl -lgslcblas -lm -o code10_gsl
//or: g++ code10_gsl.cpp -Wall -I/opt/local/include -L/opt/local/lib  -lgsl -lgslcblas -lm -o code10_gsl
//or: g++ code10_gsl.cpp -Wall -I/usr/local/include -L/usr/local/lib  -lgsl -lgslcblas -lm -o code10_gsl
//to run: ./code10_gsl


int main ()
{

  int s1 = 2;
  int s2 = 3; 

  gsl_matrix * m1 = gsl_matrix_calloc (s1, s2);

  //m1->size1
  //m2->size2
  //m1->tda //stride  实际列数
  //m1->data

  printf("tda: %ld\n",m1->tda);

  int ix =0;
  for (int i = 0; i<s1; i++)
  {
  	for (int j = 0; j<s2; j++)
  	{
      //i and j"th entry of m1
  		m1->data[i * m1->tda + j] = ix;
  		ix++;
  	}
  }


  for (int i = 0; i<s1; i++)
  {
  	for (int j = 0; j<s2; j++)
  	{
  		printf("%f ", m1->data[i * m1->tda + j]);
  	}
  	printf("\n");
  }

  //clean up
  gsl_matrix_free(m1);

  return 0;
}



