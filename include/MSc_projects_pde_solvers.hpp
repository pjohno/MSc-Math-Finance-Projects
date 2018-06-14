#pragma once
#include <vector>

namespace MSC_PROJECTS
{
// A generic lagrange interpolation function
double lagrangeInterpolation(const std::vector<double>& y,const std::vector<double>& x,double x0,unsigned int n=4);

std::vector<double> tridag(const std::vector<double>& a,const std::vector<double>& beta,const std::vector<double>& c,std::vector<double>& d);

int calculateExpectedUtility(
  double xMin, // minimum fund value
  double xMax, // minimum fund value 
  double t,  // current time t (on exit)
double T,  // final time T (on entry)
double r,  // interest rate
double mu, // risky asset returns
double sigma, // volatility of risky asset
double gamma,
int n, // grid size in fund
std::vector<double> &X,
std::vector<double> &pt,
std::vector<double> &Jt
);

}
