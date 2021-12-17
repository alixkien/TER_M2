#include <string>
#include <iostream>
#include "ressort.h"
#include "fonctions.h"

using namespace std;

Eigen::VectorXd compute_metrique(const Eigen::VectorXd & V, const Eigen::VectorXd & X, const int & n)
{
  Eigen::VectorXd m(n+2); m.setZero();
  Eigen::VectorXd u_sec(n+2);

  u_sec=derivee_seconde(V,X,n);

  for (int k=0; k<n+2; k++) {
    m(k)=max(0.5,abs(u_sec(k))); // max de u'' et de 0.5
    m(k)=sqrt(m(k));
  }
  return m;
}

Eigen::Matrix2d met_elem(const Eigen::VectorXd & M, const int & it)
{
  Eigen::Matrix2d Me; Me.setZero();
  double ki;
  ki=(M(it)+M(it+1))/2.;
  Me(0,0)+= ki;
  Me(0,1)+= -ki;
  Me(1,0)+= -ki;
  Me(1,1)+= ki;

  return Me;
}

Eigen::MatrixXd assemblage_metrique(const Eigen::VectorXd & metrique, const int & n)
{
  Eigen::MatrixXd Met(n,n); Met.setZero();
  Eigen::Matrix2d Melem; Melem.setZero();
  int I,J;
  for (int elem=0;elem<n-1;elem++) {
    Melem=met_elem(metrique,elem);
    for (int i=0; i<2; i++) {
      for (int j=0; j<2; j++) {
        I=elem+i;
        J=elem+j;
        Met(I,J)+=Melem(i,j);
      }
    }
  }
  Met(0,0)=1e12;
  Met(n-1,n-1)=1e12;

  return Met;
}
