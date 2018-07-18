#pragma once

#include <cmath>
#include <vector>
#include <functional>

namespace MSC_PROJECTS
{
    
    class QUAD
    {
        
        static constexpr double pi = M_PI;
        
        // the function A
        static double A(double x, double r, double sigma, double k, double dt);
        
        // the function B
        static double B(double x, double y, double sigma, double k, double dt);
        
        // integration proceedure: integrate the function f_i at the points y_i, contained in a vector integrand
        // PLEASE note that here we are assuming that the grid is equally spaced
        // nDown, nUp reference the subsection of the vector y we wish to integrate over
        static double integrate(int nDown, int nUp,const std::vector<double> &y,const std::vector<double> &integrand);
        
        
    public:
        
        // reset, resize and generate new grid with n+1 nodes, equally spaced on the range [min,max]
        static void resetGrid(int n, double min, double max, std::vector<double> &x);
        
        static double valueOption(double x0, double r, double sigma, double T, int n, const std::vector< double >& y, const std::vector< double >& payoff_vec)
;
        
    };
    
}
