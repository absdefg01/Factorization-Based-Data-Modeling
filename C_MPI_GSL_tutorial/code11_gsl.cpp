#include <stdio.h>
#include <gsl/gsl_matrix.h>

//to compile: g++ code11_gsl.cpp -Wall -lgsl -lgslcblas -lm -o code11_gsl
//or: g++ code11_gsl.cpp -Wall -I/opt/local/include -L/opt/local/lib  -lgsl -lgslcblas -lm -o code11_gsl
//or: g++ code11_gsl.cpp -Wall -I/usr/local/include -L/usr/local/lib  -lgsl -lgslcblas -lm -o code11_gsl
//to run: ./code11_gsl


int main ()
{

  int s1 = 2;
  int s2 = 3; 

  gsl_matrix * m1 = gsl_matrix_calloc (s1, s2);
  gsl_matrix * m2 = gsl_matrix_calloc (s1, s2);



  int ix =0;
  for (int i = 0; i<s1; i++)
  {
  	for (int j = 0; j<s2; j++)
  	{
  		m1->data[i * m1->tda + j] = ix;
      m2->data[i * m2->tda + j] = 2*ix;
  		ix++;
  	}
  }


  printf("M1\n");
  for (int i = 0; i<s1; i++)
  {
  	for (int j = 0; j<s2; j++)
  	{
  		printf("%f ", m1->data[i * m1->tda + j]);
  	}
  	printf("\n");
  }


  printf("M2\n");
  for (int i = 0; i<s1; i++)
  {
    for (int j = 0; j<s2; j++)
    {
      printf("%f ", m2->data[i * m2->tda + j]);
    }
    printf("\n");
  }

  //Scale a matrix with a constant and over-write it
  double c = 2;
  gsl_matrix_scale (m1, c); //multiply it by c

  printf("M1 Scaled\n");
  for (int i = 0; i<s1; i++)
  {
    for (int j = 0; j<s2; j++)
    {
      printf("%f ", m1->data[i * m1->tda + j]);
    }
    printf("\n");
  }
  
  //Subtract two matrices: m1 <- m1-m2 
  gsl_matrix_sub (m1, m2);

  printf("M1 Scaled - M2\n");
  for (int i = 0; i<s1; i++)
  {
    for (int j = 0; j<s2; j++)
    {
      printf("%f ", m1->data[i * m1->tda + j]);
    }
    printf("\n");
  } 


  gsl_matrix_free(m1);
  gsl_matrix_free(m2);

  return 0;
}



