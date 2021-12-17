#include <string>
#include <iostream>
#include "fonctions.h"

using namespace std;

// ========== Lecture et Sauvegarde de fichiers ========== //

Eigen::VectorXd Lecture_valeur_rho(const int & iteration, const int & nb_maille)
{
  double rho_i, temp_i, i;
  Eigen::VectorXd rho(nb_maille+2);

  //nom du fichier
  char num[10];
  sprintf(num, "%d", iteration);
  string name("solution_"), extension(".txt");
  name+=num+extension;

  // lecture du fichier
  ifstream mon_flux(name);
  for (int ligne =0; ligne <nb_maille+2; ligne ++) {
    mon_flux >> i >> rho_i >> temp_i;
    rho(ligne)=rho_i;
  }

  return rho;
}

void save_fichier(const Eigen::VectorXd maille, const Eigen::VectorXd solution, const int & nb_maille, const string & name)
{
  ofstream mon_flux; // Contruit un objet "ofstream"
  system("mkdir -p ./Resultats");

  mon_flux.open("./Resultats/" + name + ".txt");
  for (int i=0; i<nb_maille;i++) {
    mon_flux << maille(i) << " " << solution(i) << endl ;
  }

  mon_flux.close();
}

// ========== Fonctions pour le probleme ========== //

double function_Arrhenius(const double & Aref, const double & Ta, const double & rho_p, const double & T, const double & rho)
{
  return -Aref*exp(-Ta/T)*(rho-rho_p);
}

Eigen::VectorXd derivee_seconde(Eigen::VectorXd fX, Eigen::VectorXd X, const int & n)
{
  // renvoie la derivee seconde d'une fonction pour un maillage non equireparti
  Eigen::VectorXd derivee2f(n+2);
  double h1,h2;
  double c=0.5,numun=100;
  derivee2f(0)=0.0;

  for (int i=1; i<n+1; i++) {
    h1=X(i+1)-X(i);
    h2=X(i)-X(i-1);

    derivee2f(i)=2*(h2*fX(i+1)+h1*fX(i-1)-(h2+h1)*fX(i))/(h1*h1*h2+h2*h2*h1);
  }

  derivee2f(n+1)= (c*numun)*(c*numun)*(exp(1*c*numun))/(c*(1-exp(c*numun)));

  return derivee2f;
}

// ========== Resolution systeme AX=B ========== //

// Resolution LU pour une matrice A tridiagonale
Eigen::VectorXd reslu(const Eigen::MatrixXd & A, const Eigen::VectorXd & B, const int & n)
{
  Eigen::MatrixXd L(n,n), U(n,n); L.setIdentity(); U.setZero();

  // Decomposition LU sachant A tridiagonale
  U(0,0)=A(0,0); U(0,1)=A(0,1);
  for (int i=1; i<n-1; i++) {
    L(i,i-1)=A(i,i-1)/U(i-1,i-1);
    U(i,i)=A(i,i)-L(i,i-1)*A(i-1,i);
    U(i,i+1)=A(i,i+1);
  }
  L(n-1,n-2)=A(n-1,n-2)/U(n-2,n-2);
  U(n-1,n-1)=A(n-1,n-1)-L(n-1,n-2)*A(n-2,n-1);

  // Calcul de y : Ly=B
  Eigen::VectorXd vecY(n); vecY.setZero();
  vecY(0)=B(0);
  for (int k=1; k<n; k++) {
    vecY(k)=B(k)-L(k,k-1)*vecY(k-1);
  }

  // Calcul de x : Ux=y
  Eigen::VectorXd vecX(n); vecX.setZero();
  vecX(n-1)=vecY(n-1)/U(n-1,n-1);
  for (int k=n-2; k>-1; k--) {
    vecX(k)=(vecY(k)-U(k,k+1)*vecX(k+1))/U(k,k);
  }

  return vecX;
}
