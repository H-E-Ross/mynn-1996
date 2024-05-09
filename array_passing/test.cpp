i#include <stdio.h>
void cfun(double **x,const int len)
{
  printf("%d\n", len);
  printf("This is in C function cfun...\n");


  for(int i=0; i<len; i++)
    {
      printf(" %d\n  %d\n  %d\n", x[0][i]);

    }
}

void vec(  double *r[],const int len )
{
  printf("This is in C function vec...\n");

  printf("%d\n", len);

  for(int i=0; i<len; i++)
    {
       printf(" %d\n", r[i]);

    }

}
