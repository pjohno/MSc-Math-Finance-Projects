#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include "MSc_projects_pde_solvers.hpp"
#include "MSc_projects_gnuplot.hpp"
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
  calculateExpectedUtility(xMin,xMax,0.,T,r,mu,sigma,gamma,500,X,pt,Jt);

  GnuplotWidget G;
  gnuplotImage im = G.plotData(X,pt);
  cout << im.imageText << endl;
  
}

