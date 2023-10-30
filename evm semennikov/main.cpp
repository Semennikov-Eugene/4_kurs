#include <iostream>
#include "matrix.h"

int main(int argc, char const *argv[])
{
  double mu;
  double C_rho;
  double tau;
  double h_x;
  int M_x;
  int N;
  int mode;

  if (argc != 6)
    {
      printf ("Wrong number of initial values!\n");
      return -1;
    }

  if ((sscanf (argv[1], "%lf", &mu) != 1) ||
      (sscanf (argv[2], "%lf", &C_rho) != 1) ||
      (sscanf (argv[3], "%d", &N) != 1) ||
      (sscanf (argv[4], "%d", &M_x) != 1) ||
      (sscanf (argv[5], "%d", &mode) != 1)
      )
    {
      printf ("Can`t read initial values!\n");
      return -2;
    }
    h_x = 1 / M_x;
    tau = 1 / N;

    P_gas gas(1, 1, C_rho, 1.4, mu, mode);
    P_scheme scheme(M_x, N, h_x, tau);





  return 0;
}
