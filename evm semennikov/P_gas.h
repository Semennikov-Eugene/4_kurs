#pragma once
class P_gas
{
  public :
  // Параметры сетки : T и Х
  double Segm_T;
  double Segm_X;
  // Параметры газа
  double p_ro; // С_rho для зависимости C_rho*rho
  double p_gamma; // gamma для rho^gamma
  double mu; // вязкость газа mu
  bool mode; // 0 - rho^gamma; 1 - C * rho

  P_gas (double Segm_T, double Segm_X, double p_ro, double p_gamma, double mu, bool mode) :
    Segm_T(Segm_T),
    Segm_X(Segm_X),
    p_ro(p_ro),
    p_gamma(p_gamma),
    mu(mu),
    mode(mode)
  {}


};
