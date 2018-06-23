#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

// return the intesity process lambda for a given price (inverse of the price function)
double lambdaIntensity(double price,double a,double alpha)
{
  return a*exp(-alpha*price);
}

int main()
{
  // Ns is the number of stock units that may be sold
  int Ns=25;
  // parameters for the pricing function
  double a=40.,alpha=1.;
  
  // number of timesteps
  int K=100;
  // length of time T in which we can sell units
  double T=1.,dt;
  // grid variables
  dt = T/K;
  
  // create storage to track the optimal sales policy p and value function J
  // assume here that J(n,t) = J(n,t^j) = J_n^j
  // in code J[n][j] ~~ J_n^j
  vector<vector<double> > Jvalue(Ns+1,vector<double> (K+1)),optimalPvalue(Ns+1,vector<double> (K+1));
  
  // assign boundary conditions
  // J(n,T)=0
  for(int n=0;n<=Ns;n++)
  {
    Jvalue[n][K]=0.;
  }
  
  // J(0,t)=0
  for(int j=0;j<=K;j++)
  {
    Jvalue[0][j]=0.;
  }
  
  // assume that selling price is fixed p=1, we solve backwards in time so that
  //
  //  dJ/dt - lambda(p) ( J(n,t+dt)-J(n-1,t+dt) ) + p lambda(p) = 0
  // and in finite differences we have
  //  ( J(n,t+dt) - J(n,t) )/dt - lambda(p) ( J(n,t+dt)-J(n-1,t+dt) ) + p lambda(p) = 0
  // then we can write
  //      J(n,t) = ( 1 - lambda(p) dt ) J(n,t+dt) + ( lambda(p) dt ) J(n-1,t+dt) + p lambda(p)
  //      J_n^j = ( 1 - lambda(p) dt ) J_n^{j+1} + ( lambda(p) dt ) J_{n-1}^{j+1} + p lambda(p)
  // notice the similarity to binomial tree pricing formula, except q=lambda(p) dt is the probability
  //
  // At each time step, calculate the Jvalue
  for(int j=K-1;j>=0;j--)
  {
    for(int n=1;n<=Ns;n++)
    {
      double pStar = Jvalue[n][j+1] - Jvalue[n-1][j+1] + alpha;
      double lambdaStar = lambdaIntensity(pStar,a,alpha);
      Jvalue[n][j]=(1.- lambdaStar*dt) * Jvalue[n][j+1] + 
	lambdaStar*dt * Jvalue[n-1][j+1] +
	pStar * lambdaStar * dt;
      optimalPvalue[n][j] = pStar;
    }
  }
  // output values to screen
  for(int n=0;n<=Ns;n++)
  {
    for(int j=0;j<=K;j++)
    {
      cout << n << " " << j*dt << " " << Jvalue[n][j] << " " << optimalPvalue[n][j] <<  endl;
    }
    cout << endl;
  }
}
