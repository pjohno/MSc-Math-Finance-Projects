#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>

const double pi = 4.*atan(1.);
// the function A
double  A(double x, double r, double sigma, double k, double dt)
{
    return 1. / sqrt(2 * sigma*sigma*pi*dt) * exp(-0.5*k*x - 0.125*sigma*sigma*k*k*dt - r*dt);
}

// the function B
double  B(double x, double y, double sigma, double k, double dt)
{
    return exp(-(x - y)*(x - y) / (2.*sigma*sigma*dt) + 0.5 * k * y);
}

// integration proceedure: integrate the function f_i at the points y_i, contained in a vector integrand
// PLEASE note that here we are assuming that the grid is equally spaced
// nDown, nUp reference the subsection of the vector y we wish to integrate over
double  integrate(int nDown, int nUp,const std::vector<double> &y,const std::vector<double> &integrand)
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
void  resetGrid(int n, double min, double max, std::vector<double> &x)
{
    x.resize(n + 1);
    double dx = (max - min) / double(n);
    if (n == 0)dx = 0.;
    for (uint i = 0; i<x.size(); i++) x[i] = min + i*dx;
}

double  valueOption(double  S0,double r,double sigma,double T,int n,const std::vector<double> &y,const std::vector<double> &payoff_vec)
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

double callOption(double  S,double X,double r,double sigma,double T,int N=500,double xi=7.5)
{
    // set y grid up at the final step, using strike price as lower limit
    std::vector<double> y(N+1),V(N+1);
    resetGrid(N,log(X),log(X*exp(xi*sigma*sqrt(T))),y);
    // assign payoff, since e^y>=X we don't need max function
    for(unsigned int i=0;i<y.size();i++)
        V[i] = exp(y[i])-X;
    
    // return the value V(x=log(S),0) = A(x)*\int B(x,y)*V(y,T) dy
    return valueOption(S,r,sigma,T,N,y,V);
}

double putOption(double  S,double X,double r,double sigma,double T,int N=500,double xi=7.5)
{
    // set y grid up at the final step, using strike price as lower limit
    std::vector<double> y(N+1),V(N+1);
    double yMin = log(X*exp(-xi*sigma*sqrt(T)));
    double yMax = log(X);
    double dy = (yMax - yMin)/N;
    for(int i=0;i<=N;i++)
        y[i] = yMin + i*dy;
    // assign payoff, since e^y>=X we don't need max function
    for(unsigned int i=0;i<y.size();i++)
        V[i] = X-exp(y[i]);
    
    // return the value V(x=log(S),0) = A(x)*\int B(x,y)*V(y,T) dy
    return valueOption(S,r,sigma,T,N,y,V);
}

int main()
{
    double  S = 100., X = 100., r = 0.06, sigma = 0.2, T = 1.;
    std::cout << " A call option using QUAD...\n\n";
    std::cout << "N     |   V(S,0;N) |   R   |   conv rate" << std::endl;
    double valueOld = 1.,diffOld=1.;
    int k=2;
    for(int n=10;n<=5000;n*=k)
    {
        double value = putOption(S,X,r,sigma,T,n);
        double diff = value - valueOld;
        double R = diffOld/diff;
        std::cout <<  n << " | " << value << " | " << R << " | " << log(R)/log(k) << std::endl;
        valueOld = value;
        diffOld = diff;
    }
    
}
