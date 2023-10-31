#include "matrix.h"

double d_rho (double x, double t)
{
  return -3 * M_PI * exp (t) * sin (3 * x * M_PI);
}

double Matrix:: get_Fn_m (int n, int m, bool mode)
{
  double x = scheme.h_x * m;
  double t = scheme.tau * n;

  double result = (-2 * M_PI * sin (2 * M_PI * t)) * sin (4 * M_PI * x) * rho (x, t)
                  + rho (x, t) * u (x, t) * cos (2 * M_PI * t) * 4 * M_PI * cos (4 * M_PI * x)
                  + gas.mu * 16 * M_PI * M_PI * u (x,t) + get_Pn_m(n, m, mode);
  return result;
}

double Matrix:: get_Pn_m (int n, int m, bool mode)
{
  double x = scheme.h_x * m;
  double t = scheme.tau * n;

  double result = mode ? gas.p_ro * d_rho(x, t) :
                  gas.p_gamma * pow( rho(x, t), gas.p_gamma - 1) * d_rho(x, t);
  return result;
}


void Matrix::init_vector_G ()
{
  double tau = scheme.tau;
  double h_x = scheme.h_x;

  printf("solution_V\n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", solution_V[i]);

  printf("\nsolution_G\n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", solution_G[i]);

  // основная часть
  for (int m = 1; m < Dim-1; m++)
    {
      double Vn_m = solution_V[m];
      double Vn_m_r = solution_V[m+1];
      double Vn_m_l = solution_V[m-1];
      matrix_vec_a[m] = 1;
      matrix_vec_b[m] = (tau / (4 * h_x)) * (Vn_m + Vn_m_r);
      matrix_vec_c[m] = -(tau / (4 * h_x)) * (Vn_m + Vn_m_l);

      double Gn_m = solution_G[m];
      rhs_vector[m] = Gn_m + (tau * (Gn_m - 2) * (Vn_m_r - Vn_m_l)) / (2 * h_x);
    }

  //граничные условия : 0 строка
  double Vn_0 = solution_V[0];
  double Vn_1 = solution_V[1];
  double Vn_2 = solution_V[2];
  double Vn_3 = solution_V[3];
  double Gn_0 = solution_G[0];
  double Gn_1 = solution_G[1];
  double Gn_2 = solution_G[2];
  double Gn_3 = solution_G[3];

  matrix_vec_c[0] = 0;
  matrix_vec_a[0] = (1 / tau - Vn_0 / (2 * h_x));
  matrix_vec_b[0] = Vn_1 / (2 * h_x);
  rhs_vector[0] = Gn_0 / tau + 0.5 * ((Gn_0 - 2) * (Vn_1 - Vn_0) / (h_x)) +
                  0.5 * ((4 * Gn_2 * Vn_2 - 5 * Gn_1 * Vn_1 + 2 * Gn_0 * Vn_0 - Gn_3 * Vn_3 ) / (2 * h_x) +
                  (2 - Gn_0) * (4 * Vn_2 - 5 * Vn_1 + 2 * Vn_0 - Vn_3) / (2 * h_x));

  // граничные условия : строка Dim-1
  double Vn_M_0 = solution_V[Dim - 1];
  double Vn_M_1 = solution_V[Dim - 2];
  double Vn_M_2 = solution_V[Dim - 3];
  double Vn_M_3 = solution_V[Dim - 4];
  double Gn_M_0 = solution_G[Dim - 1];
  double Gn_M_1 = solution_G[Dim - 2];
  double Gn_M_2 = solution_G[Dim - 3];
  double Gn_M_3 = solution_G[Dim - 4];

  matrix_vec_b[Dim - 1] = 0;
  matrix_vec_a[Dim - 1] = (1 / tau - Vn_M_0 / (2 * h_x));
  matrix_vec_c[Dim - 1] = -Vn_M_1 / (2 * h_x);
  rhs_vector[Dim - 1] = Gn_M_0 / tau + 0.5 * ((Gn_M_0 - 2) * (Vn_M_0 - Vn_M_1) / (h_x)) -
                  0.5 * ((4 * Gn_M_2 * Vn_M_2 - 5 * Gn_M_1 * Vn_M_1 + 2 * Gn_M_0 * Vn_M_0 - Gn_M_3 * Vn_M_3 ) / (2 * h_x) +
                  (2 - Gn_M_0) * (4 * Vn_M_2 - 5 * Vn_M_1 + 2 * Vn_M_0 - Vn_M_3) / (2 * h_x));

}

double Matrix::get_mu ()
{
  double max = -1;

  for (int i = 0; i < Dim; i++)
    {
      double candidate = exp(-solution_G[i]);
      if (candidate > max)
        max = candidate;
    }
  return gas.mu * max;
}

void Matrix::init_vector_V ()
{
  // первая строка фиктивная
  matrix_vec_a[0] = 1;
  // printf("AAAAAAAaa%lf\n", matrix_vec_a[0]);
  matrix_vec_b[0] = 0;
  matrix_vec_c[0] = 0;
  rhs_vector[0] = 0;

  // основная часть
  double tau = scheme.tau;
  double h_x = scheme.h_x;
  double mu = get_mu();
  for (int m = 1; m < Dim-1; m++)
    {
      matrix_vec_a[m] = (1 / tau + (2 * mu) / (h_x * h_x));

      double Vn_m = solution_V[m];
      double Vn_m_r = solution_V[m + 1];
      double Vn_m_l = solution_V[m-1];

      matrix_vec_b[m] = (Vn_m + Vn_m_r) / (6 * h_x) - mu / (h_x * h_x);
      matrix_vec_c[m] = -(Vn_m + Vn_m_l) / (6 * h_x) - mu / (h_x * h_x);

      double Gn_m = solution_G[m];
      double Gn_m_r = solution_G[m + 1];
      double Gn_m_l = solution_G[m - 1];
      double Pn_m = get_Pn_m(step, m, gas.mode);
      double Fn_m = get_Fn_m(step, m, gas.mode);

      rhs_vector[m] = Vn_m / tau  - Pn_m * Gn_m * (Gn_m_r - Gn_m_l) / (2 * h_x)
                      - (mu - gas.mu * exp (Gn_m)) * (Vn_m_r - 2 * Vn_m + Vn_m_l)
                       / (h_x * h_x);
    }

    // printf ("")

  // последняя строка тоже фиктивная
  matrix_vec_a[Dim - 1] = 1;
  matrix_vec_b[Dim - 1] = 0;
  matrix_vec_c[Dim - 1] = 0;
  rhs_vector[Dim - 1] = 0;
}


void Matrix::solve (bool mode)
{
  std::vector<double> & solution = mode ? solution_V : solution_G;
  // Метод прогонки : вперед

  printf ("\nsolving matrix mdoe %d\n", mode);
  printf ("\n matrix a \n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", matrix_vec_a[i]);

  printf ("\n\n\n\n matrix b \n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", matrix_vec_b[i]);

  printf ("\n\n\n\n matrix c \n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", matrix_vec_c[i]);

  printf ("\n\n\n\n rhs \n");
  for (int i = 0; i< Dim; i++)
    printf("%lf ", rhs_vector[i]);

  printf ("\n\n");

  double znam = matrix_vec_a[0];
  matrix_vec_a[0] =  (-1 * matrix_vec_b[0]) / znam;
  matrix_vec_b[0] = rhs_vector[0] / znam;

  for (int i = 1; i < Dim - 1; i++)
    {
      znam = matrix_vec_a[i] + matrix_vec_c[i] * matrix_vec_a[i-1];
      matrix_vec_a[i] = (-1 * matrix_vec_b[i]) / znam;
      double chis = rhs_vector[i] - matrix_vec_c[i] * matrix_vec_b[i-1];
      matrix_vec_b[i] = chis / znam;
    }
  // Метод прогонки : обратно
  solution[Dim - 1] = (rhs_vector[Dim - 1] - matrix_vec_c[Dim - 1]*matrix_vec_b[Dim - 2]) /
                      (matrix_vec_a[Dim - 1] + matrix_vec_c[Dim - 1] * matrix_vec_a[Dim - 2]);
  for (int i = Dim - 2; i >=0; i--)
    {
      solution[i] = matrix_vec_a[i] * solution[i+1] + matrix_vec_b[i];
    }

   printf ("\n\nsolution mode = %d\n", mode);
    for (int i = 0; i< Dim; i++)
      printf ("%lf ", solution [i]);
    printf ("\n\n");
}
