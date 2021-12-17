#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen"

Eigen::VectorXd Lecture_valeur_rho(const int & iteration, const int & nb_maille);

void save_fichier(const Eigen::VectorXd maille, const Eigen::VectorXd solution, const int & nb_maille, const std::string & name);

Eigen::VectorXd derivee_seconde(Eigen::VectorXd fX, Eigen::VectorXd X, const int & n);

double function_Arrhenius(const double & Aref, const double & Ta, const double & rho_p, const double & T, const double & rho);

Eigen::VectorXd reslu(const Eigen::MatrixXd & A, const Eigen::VectorXd & B, const int & n);
