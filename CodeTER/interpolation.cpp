#include <string>
#include <iostream>
#include "interpolation.h"
#include "fonctions.h"

using namespace Eigen;
using namespace std;

VectorXd resol_sys_temp_adapt(VectorXd X, VectorXd oldX, VectorXd T, VectorXd lambda, VectorXd rho, VectorXd Cp, double dt, double Lm, double A_ref, double rho_p, double Ta)
{
  int N=T.rows();
  MatrixXd A(N,N);              // matrice du système
  VectorXd hn1(N);              //h n+1
  VectorXd hn(N);               //h n
  VectorXd f(N);                // second membre f tilde
  VectorXd gamma(N);            //vecteur gamma
  double phi=1e5;               //flux
  double alphai;
  double nu1;                   //vi+1/2
  double nu2;                   //vi-1/2
  A.setZero(N,N);

  //calcul hn et hn+1
  for (int i=0; i<N; i++) {
    hn1(i)=(X(i+2)-X(i))/2;
    hn(i)=(oldX(i+2)-oldX(i))/2;
  }

  //premiere ligne de A
  A(0,0)=1 + 2*lambda(0)*dt/(rho(0)*Cp(0)*hn1(0)*(hn1(0)+hn1(1)));
  A(0,1)= -2*lambda(0)*dt/(rho(0)*Cp(0)*hn1(0)*(hn1(0)+hn1(1)));

  //calcul ligne 1 de f tilde
  gamma(0)= (dt*Lm)/(rho(0)*Cp(0))*(-function_Arrhenius(A_ref, Ta,rho_p,T(0),rho(0)));
  nu1=((X(2)+oldX(2)) - (X(1)+oldX(1)))/(dt);
  nu2=((X(1)+oldX(1)) - (X(0)+oldX(0)))/(dt);
  f(0)=gamma(0)+hn1(0)/hn(0)*T(0)+(dt/hn1(0))*(fmax(nu1,0)*T(0)+fmin(nu1,0)*T(1)-fmax(nu2,0)*T(0)-fmin(nu1,0)*(T(0)+(phi*hn1(0)/lambda(0))));

  for (int i=1; i<N-1; i++)
  {
    //calcul de A
    alphai=2*lambda(i)/((rho(i)*Cp(i)*hn1(i))*(hn1(i)+hn1(i+1))*(hn1(i)+hn1(i-1)));
    A(i,i-1) = -alphai*(hn1(i+1)+hn1(i));
    A(i,i)   = 1 + alphai*(hn1(i-1)+2*hn1(i)+hn1(i+1));
    A(i,i+1) = -alphai*(hn1(i)+hn1(i-1));

    //calcul de f tilde
    gamma(i)= Lm*dt/(rho(i)*Cp(i))*(-function_Arrhenius(A_ref, Ta,rho_p,T(i),rho(i)));
    nu1=((X(i+2)+oldX(i+2)) - (X(i+1)+oldX(i+1)))/(dt);
    nu2=((X(i+1)+oldX(i+1)) - (X(i)+oldX(i)))/(dt);
    f(i)=gamma(i)+hn1(i)/hn(i)*T(i)+(dt/hn1(i))*(fmax(nu1,0)*T(i)+fmin(nu1,0)*T(i+1)-fmax(nu2,0)*T(i-1)-fmin(nu1,0)*T(i));

  }

  //derniere ligne de A
  A(N-1,N-2)= -2*lambda(N-1)*dt/(rho(N-1)*Cp(N-1)*hn1(N-1)*(hn1(N-1)+hn1(N-2)));
  A(N-1,N-1)= 1 + 2*lambda(N-1)*dt/(rho(N-1)*Cp(N-1)*hn1(N-1)*(hn1(N-1)+hn1(N-2)));

  //derniere ligne de f tilde
  gamma(N-1)=Lm*dt/(rho(N-1)*Cp(N-1))*(-function_Arrhenius(A_ref, Ta,rho_p,T(N-1),rho(N-1)));
  nu1=((X(N+1)+oldX(N+1)) - (X(N)+oldX(N)))/(dt); //vi +1/2
  nu2=((X(N)+oldX(N)) - (X(N-1)+oldX(N-1)))/(dt); //vi -1/2
  f(N-1)=gamma(N-1)+hn1(N-1)/hn(N-1)*T(N-1)+(dt/hn1(N-1))*(fmax(nu1,0)*T(N-1)+fmin(nu1,0)*T(N-1)-fmax(nu2,0)*T(N-2)-fmin(nu1,0)*T(N-1));

  // resout AT(n+1)=f
  T=reslu(A,f,N);

  return T;
}

VectorXd resol_sys_temp(VectorXd X, VectorXd T, VectorXd lambda, VectorXd rho, VectorXd Cp, double dt, double Lm, double A_ref, double rho_p, double Ta)
{
  int N=T.rows();      //N=taille de T
  MatrixXd A(N,N);     // matrice du système
  VectorXd h(N);       //h
  VectorXd f(N);       // second membre
  double phi=1e5;      //flux
  double alphai;
  A.setZero(N,N);

  //calcul de h
  for (int i=0; i<N; i++)
  {
    h(i)=(X(i+2)-X(i))/2;
  }

  //premiere ligne de A
  A(0,0)=1 + 2*lambda(0)*dt/(rho(0)*Cp(0)*h(0)*(h(0)+h(1)));
  A(0,1)= -2*lambda(0)*dt/(rho(0)*Cp(0)*h(0)*(h(0)+h(1)));

  //premiere ligne de f
  f(0)= (dt*Lm)/(rho(0)*Cp(0))*(-function_Arrhenius(A_ref, Ta,rho_p,T(0),rho(0))) + (dt*phi)/(rho(0)*Cp(0)*h(0));

  for (int i=1; i<N-1; i++)
   {
    //matrice A
    alphai=2*lambda(i)/((rho(i)*Cp(i)*h(i))*(h(i)+h(i+1))*(h(i)+h(i-1)));
    A(i,i-1) = -alphai*(h(i+1)+h(i));
    A(i,i)   = 1 + alphai*(h(i-1)+2*h(i)+h(i+1));
    A(i,i+1) = -alphai*(h(i)+h(i-1));

    //second membre f
    f(i)= Lm*dt/(rho(i)*Cp(i))*(-function_Arrhenius(A_ref, Ta,rho_p,T(i),rho(i)));
   }

   //derniere ligne de A
  A(N-1,N-2)= -2*lambda(N-1)*dt/(rho(N-1)*Cp(N-1)*h(N-1)*(h(N-1)+h(N-2)));
  A(N-1,N-1)= 1 + 2*lambda(N-1)*dt/(rho(N-1)*Cp(N-1)*h(N-1)*(h(N-1)+h(N-2)));

  //derniere ligne second membre
  f(N-1)=Lm*dt/(rho(N-1)*Cp(N-1))*(-function_Arrhenius(A_ref, Ta,rho_p,T(N-1),rho(N-1)));

  // resout AT(n+1)=(f+T(n))
  T=reslu(A,T+f,N);

  return T;
}
