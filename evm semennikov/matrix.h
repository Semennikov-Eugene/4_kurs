#pragma once
#include <cmath>
#include "P_scheme.h"
#include "P_gas.h"
#include <vector>

enum class mode : bool {V = true, G = false};

class Matrix
{
  double rho_0 (double x)
  {
    return cos(3 * M_PI * x) + 1.5;
  }

  double u_0 (double x)
  {
    return sin(4 * M_PI * x);
  }

  double get_Pn_m (int n, int m, bool mode);
  double get_Fn_m (int n, int m, bool mode);

  P_gas gas;
  P_scheme scheme;

  int Dim;
  int step = 0;
  std::vector<double> matrix_vec_a;
  std::vector<double> matrix_vec_b;
  std::vector<double> matrix_vec_c;
  std::vector<double> rhs_vector;
  std::vector<double> solution_G;
  std::vector<double> solution_V;

  Matrix (P_gas gas, P_scheme scheme) :
  gas(gas), scheme(scheme)
  {
    Dim = scheme.Dim;
    matrix_vec_a.resize(Dim);
    matrix_vec_b.resize(Dim);
    matrix_vec_c.resize(Dim);
    rhs_vector.resize(Dim);
    solution_G.resize(Dim);
    solution_V.resize(Dim);
    for (int i = 0; i < Dim; i++)
      {
        double x = scheme.get_point_x_by_i(i);
        solution_G[i] = log(rho_0(x));
        solution_V[i] = u_0(x);
      }
    step = 0;
  }

  double get_mu ();
  void init_vector_G();
  void init_vector_V();
  void solve (bool mode /*0 if solution for G 1 otherwise*/);
};
