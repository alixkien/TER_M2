#include <string>
#include "fonctions.h"
#include "adapt.h"
#include "interpolation.h"

using namespace Eigen;
using namespace std;

//Appliquer fonction "fonc" pour avoir la solution rho aux points du maillage
//calculer rho'' en appliquant la fonction rho_seconde avec le vecteur rho (déja calculé aux différents points)
//on applique le changement de maillage avec la fonction adaptation
//on applique interpole pour avoir les valeurs de rho aux nouveaux points
//etc

int main()
{
  int N(100), iteration_max(101);
  VectorXd X(N+2), oldX(N+2), newX(N), rho(N), rho_etoile(N), T(N);
  X.setZero(); rho.setZero();

  double L=1.;  // longueur du maillage
  double dx=L/(N+1); // pas d'espace
  double dt=0.1; // pas de temps

  // definition des constantes physiques du probleme
  double Aref(1000.),Ta(6000.),rho_v(1500.),rho_p(1000.),Cpv(1000.),Cpp(1500.),Lm(3e6);
  //double lambda_v(1.5),lambda_p(1);

  // =====  maillage intial  ===== //
  for (int i=0; i<N+2; i++) {
    X(i)=i*dx;
  }

  // ====  initialisation de rho et de la temperature  ===== //
  for (int i=0; i<N ; i++) {
    T(i)=293; rho(i)=1500;
  }

  //====  Resolution des EDO  ===== //
  for (int iteration=0; iteration<iteration_max; iteration++) {

    cout << "//------------ iteration " << to_string(iteration) << "--------- //" << endl;

    // Calcul du nouveau maillage
    newX=adaptation(X,rho,N,L);
    if (iteration<10) {
      save_fichier(X, newX, N+2, "adapt" + to_string(iteration));
    }

    //  Etape 1: calcul de rho*
    cout << "===== Etape 1 =====" << endl;
    for (int i=0 ;i<N ;i++) {
      rho_etoile(i)=(rho(i) +dt*rho_p*Aref*exp(-Ta/T(i)))/(1+dt*Aref*exp(-Ta/T(i)));
    }

    // Etape 2-3: évaluation de T avec rho*
    cout << "===== Etape 2 =====" << endl;
    VectorXd xi(N),Cp(N),lambda(N);
    for (int i=0; i<N ; i++) {
      xi(i)=(rho_v-rho_etoile(i))/(rho_v-rho_p);
      Cp(i)=((1-xi(i))*rho_v*Cpv+xi(i)*rho_p*Cpp)/rho_etoile(i);
      //lambda(i)=(1-xi(i))*lambda_v+xi(i)*lambda_p;
      lambda(i)=1;
    }

    cout << "===== Etape 3 =====" << endl;
    T=resol_sys_temp_adapt(newX,X,T,lambda,rho_etoile,Cp,dt,Lm,Aref,rho_p,Ta);

    // Etape 4 : réévaluation de rho
    cout << "===== Etape 4 =====" << endl;
    for (int i=0 ; i<N ; i++) {
      rho(i)=rho(i)+dt*function_Arrhenius(Aref,Ta,rho_p,T(i),rho(i));
    }

    X=newX;

    cout << " - - - - - End iteration - - - - - " << endl;

    if (iteration%10==0) {
      cout << " - - - - - Save - - - - - " << endl;
      save_fichier(newX, rho , N, "rho" + to_string(iteration));
      save_fichier(newX, T , N, "temp" + to_string(iteration));
    }
  }

  return 0;
}
