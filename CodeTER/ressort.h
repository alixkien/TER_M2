#include <iostream>
#include <cmath>
#include <algorithm>    // std::max
#include "Eigen"

Eigen::VectorXd compute_metrique(const Eigen::VectorXd & V, const Eigen::VectorXd & X, const int & n);

Eigen::Matrix2d met_elem(const Eigen::VectorXd & M, const int & it);

Eigen::MatrixXd assemblage_metrique(const Eigen::VectorXd & m, const int & n);
