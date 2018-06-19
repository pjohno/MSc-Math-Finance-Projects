#include <MSc_projects_quad.hpp>
#include <algorithm>

namespace MSC_PROJECTS
{
    // the function A
    double QUAD::A(double x, double r, double sigma, double k, double dt)
    {
        return 1. / sqrt(2 * sigma*sigma*pi*dt) * exp(-0.5*k*x - 0.125*sigma*sigma*k*k*dt - r*dt);
    }
    
    // the function B
    double QUAD::B(double x, double y, double sigma, double k, double dt)
    {
        return exp(-(x - y)*(x - y) / (2.*sigma*sigma*dt) + 0.5 * k * y);
    }
    
    // integration proceedure: integrate the function f_i at the points y_i, contained in a vector integrand
    // PLEASE note that here we are assuming that the grid is equally spaced
    // nDown, nUp reference the subsection of the vector y we wish to integrate over
    double QUAD::integrate(int nDown, int nUp,const std::vector<double> &y,const std::vector<double> &integrand)
    {
        double sum, h = (y[nUp] - y[nDown]) / (nUp - nDown);
        sum = integrand[nDown];
        for (int i = nDown+2; i<nUp; i+=2)
        {
            sum += 2.*integrand[i];
        }
        for (int i = nDown+1; i<nUp; i+=2)
        {
            sum +=  4.*integrand[i];
        }
        sum += integrand[nUp];
        return h / 3.*sum;
    }
    
    // reset, resize and generate new grid with n+1 nodes, equally spaced on the range [min,max]
    void QUAD::resetGrid(int n, double min, double max, std::vector<double> &x)
    {
        x.resize(n + 1);
        double dx = (max - min) / double(n);
        if (n == 0)dx = 0.;
        for (uint i = 0; i<x.size(); i++) x[i] = min + i*dx;
    }
    
    double QUAD::valueOption(double  S0,double r,double sigma,double T,int n,const std::vector<double> &y,const std::vector<double> &payoff_vec)
    {
        // value of x
        double x = log(S0);
        
        // variable k
        double k = 2.*r / sigma / sigma - 1.;
        
        // setup vectors to store option values, and integrand function
        // these vectors must be the same size as the y vector
        std::vector<double> integrand(y.size());
        
        for (uint j = 0; j < y.size(); j++)
        {
            integrand[j] = B( x , y[j] , sigma , k , T ) * payoff_vec[j];
        }
        
        // in this example we integrate over the entire range of y, it is possible to select a subrange
        return A(x, r, sigma, k, T )*integrate(0, n, y, integrand);
        
    }
}
