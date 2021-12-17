#include <iostream>
#include <cmath>
#include <algorithm>    // std::max
#include "Eigen"

Eigen::VectorXd adaptation(const Eigen::VectorXd & maille, const Eigen::VectorXd & u, const int & n, const int & L);
