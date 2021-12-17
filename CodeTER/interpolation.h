#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen"

Eigen::VectorXd resol_sys_temp_adapt(Eigen::VectorXd X, Eigen::VectorXd oldX, Eigen::VectorXd T, Eigen::VectorXd lambda, Eigen::VectorXd rho, Eigen::VectorXd Cp, double dt, double Lm, double A_ref, double rho_p, double Ta);

Eigen::VectorXd resol_sys_temp(Eigen::VectorXd X, Eigen::VectorXd T, Eigen::VectorXd lambda, Eigen::VectorXd rho, Eigen::VectorXd Cp, double dt, double Lm, double A_ref, double rho_p, double Ta);

//Eigen::VectorXd interpole(Eigen::VectorXd xn, Eigen::VectorXd yn, Eigen::VectorXd xn1, const int & n);

//Eigen::VectorXd fonc(Eigen::VectorXd xn,Eigen::VectorXd rho, Eigen::VectorXd T);

//Eigen::VectorXd Euler_explicite(Eigen::VectorXd xn, Eigen::VectorXd yn_interpole, const double & t, const double & dt);

//Eigen::VectorXd RK4(Eigen::VectorXd xn, Eigen::VectorXd yn_interpole, const double & t, const double & dt);
