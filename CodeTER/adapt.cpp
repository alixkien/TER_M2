#include <string>
#include <iostream>
#include "adapt.h"
#include "ressort.h"
#include "fonctions.h"

using namespace std;

Eigen::VectorXd adaptation(const Eigen::VectorXd & maille, const Eigen::VectorXd & u, const int & n, const int & L)
{ // X vecteurs des points du maillage, u vecteur des solutions
  int ite(0), itemax(10);  // ite boucle pour s'assurer que ca s'arrete
  double tol(0.0001); // boucle de precision
  bool test=true;
  Eigen::MatrixXd mat_metric(n+2,n+2); // matrice de la metrique
  Eigen::VectorXd metric(n+2); // vecteur des metriques
  Eigen::VectorXd Xadapt(n+2); // Xadapt nouveau point du maillage
  Eigen::VectorXd D(n+2), B(n+2), sol(n+2), X(n+2); X=maille;
  sol(0)=u(0);
  for (int i=1; i<n+1; i++) {
    sol(i)=u(i-1);
  }
  sol(n+1)=u(n-1);

  while ((test)&&(ite<itemax)) {
    // ===== Remplissage du vecteur des metriques ===== //
    metric=compute_metrique(sol,X,n);

    // ===== Creation de la matrice de la metrique ===== //
    mat_metric=assemblage_metrique(metric,n+2);

    // ===== Definition du nouveau maillage ===== //
    B.setZero(); B(n+1)=1e12;

    // ===== Résolution A*Xadapt = B ===== //
    Xadapt=reslu(mat_metric,B,n+2);
    D.setZero(); D=mat_metric*Xadapt-B;

    D.setZero(); D=X-Xadapt;
    if (D.norm()<tol) {
      test=false;
    }
    else {
      X=Xadapt; // Mise à jour des points du maillage
      ite++;
    }
  }
  return X;
}
