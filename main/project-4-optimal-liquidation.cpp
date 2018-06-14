#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "MSc_projects_mVector.hpp"
#include "MSc_projects_Integrate.hpp"
using namespace std;
using namespace MSC_PROJECTS;

int main()
{
  
  // parameters from Gueant et al 2012
  // figure 1
  int qMax=6;
  double mu=0.0;
  double sigma=0.3;
  double A = 0.1;
  double k=0.3;
  double gamma = 0.05;
  double b=3;
  double T= 300;
  
  // number of time observations
  int n=100;
  // time step
  double dT = T/n;
  
  // store the value of the omega and delta at each time step and q value
  vector<MVector> omega(n+1,MVector(qMax+1)),delta(n+1,MVector(qMax+1));
  
  // initialise the solution at t=T
  for(int q=0;q<=qMax;q++)
  {
    // from initial condition
    omega[n][q] = exp(-k*q*b);
    // from formula for delta
    if(q>0)
      delta[n][q] = 1/k*log(omega[n][q]/omega[n][q-1]) + 1/gamma*log(1+gamma/k);
    else
      delta[n][q] = 0.;
  }
  
  // now use a numerical integration to find the value at T_i given the value at T_{i+1}
  for(int i = n-1 ; i>=0 ; i--)
  {
    // some constants in the equation
    double alpha = k/2.*gamma*sigma*sigma;
    double beta = k*mu;
    double eta = A*pow(1+gamma/k,-(1+k/gamma));
    // now solve 
    //   dw(q,t)/dt = (alpha q^2 - beta q) w(q,t) - eta w(q-1,t)
    // with initial condition 
    //   w(q,T_i) = omega(q,T_i)
    // so that
    //   omega(q,T_{i-1}) = w(q,T_{i-1})
    omega[i] = RK4MethodTemplate(100,i*dT,(i-1)*dT,omega[i+1],
				 [&]
				 (const  MVector &w,double t)
				 {
				   MVector F(qMax+1);
				   F[0] = 0.;
				   for(int q=1;q<=qMax;q++)
				     F[q] = (alpha*q*q - beta*q)*w[q] - eta*w[q-1];
				   return F;
				 }
    );
    
    // We can then calculate the value of the optimal ask price
    for(int q=1;q<=qMax;q++)
      delta[i][q] = 1/k*log(omega[i][q]/omega[i][q-1]) + 1/gamma*log(1+gamma/k);
    delta[i][0] = 0.;
  }
  
  // and output results to file
  ofstream output("test.csv");
  for(int i = 0 ; i<=n ; i++)
  {
    output << i*dT << " , " << omega[i] << " , " << delta[i] << endl;
  }
  
}

