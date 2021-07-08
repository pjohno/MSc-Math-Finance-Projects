#pragma once
#include <vector>

namespace MSC_PROJECTS
{
    
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
