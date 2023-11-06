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

    h_x = 1. / M_x;
    tau = 1. / N;

    P_gas gas (1, 1, C_rho, 1.4, mu, mode);
    P_scheme scheme (M_x, N, h_x, tau);
    Matrix matrix (gas, scheme);


//    info for vec_V will be calculated there
    matrix.init_vector_G ();
    //now solver works properly
    matrix.solve (0);
    matrix.calc_and_print_nev_C (0, 1);
    // init and solve
    matrix.init_vector_V ();
    matrix.solve (1);
    matrix.calc_and_print_nev_C (1, 1);
    // step 1 passed go futher
    for (int i = 2; i <= N; i++)
      {
        matrix.init_vector_G ();
        matrix.solve (0);
        matrix.calc_and_print_nev_C (0, i);
        matrix.init_vector_V ();
        matrix.solve (1);
        matrix.calc_and_print_nev_C (1, i);
      }

  return 0;
}
