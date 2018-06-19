#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "MSc_projects_pde_solvers.hpp"
using namespace std;
using namespace MSC_PROJECTS;

int main()
{
  // problem parameters
  double gamma=-3;
  double sigma=0.05;
  double r=0.0538;
  double mu=r+0.02;
  double T=1;
  double xMin=0.5; // minimum fund value
  double xMax=5; // minimum fund value 
  vector<double> X,pt,Jt;
  calculateExpectedUtility(xMin,xMax,0.,T,r,mu,sigma,gamma,1000,X,pt,Jt);

  std::vector<double> investmentYield(Jt.size());
  for(int i=0;i<Jt.size();i++)
  {
      // investment yield = 1/T * log(certainty equivalent / initial investment)
      investmentYield[i] = 1./T * log(pow( gamma * Jt[i] , 1./gamma ) / X[i]);
  }
  
  // and output results to file
  ofstream output("test.csv");
  for(int i = 0 ; i<X.size() ; i++)
  {
    output << X[i] << " , " << Jt[i] << " , " << pt[i]  << " , " << investmentYield[i] << endl;
  }
  
}

