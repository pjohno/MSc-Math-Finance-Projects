#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <functional>

const double m_pi = 4.*atan(1.);

// the function A
double A(double x, double r, double sigma, double k, double dt)
{
  return 1. / sqrt(2 * sigma*sigma*m_pi*dt) * exp(-0.5*k*x - 0.125*sigma*sigma*k*k*dt - r*dt);
}

// the function B
double B(double x, double y, double sigma, double k, double dt)
{
  return exp(-(x - y)*(x - y) / (2.*sigma*sigma*dt) + 0.5 * k * y);
}

// integration proceedure: integrate the function f_i at the points y_i, contained in a vector integrand
// PLEASE note that here we are assuming that the grid is equally spaced
// nDown, nUp reference the subsection of the vector y we wish to integrate over
double integrate(int nDown, int nUp,const std::vector<double> &y,const std::vector<double> &integrand)
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

double valueOption(double  x,double r,double sigma,double T,int n,
		   const std::vector<double> &y,const std::vector<double> &payoff_vec)
{
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

int main()
{
  double  S,X,r,sigma,dT;
  
  S = 76.; X = 100.; r = 0.06; sigma = 0.2; dT = 0.1;
  
  {
    int N=500;
    double x0 = log(S);
    
    std::cout << " A perpetual Bermudan put option ...\n";
    
    // for a perpetual American put option we have 
    double alpha = -2*r/sigma/sigma;
    double theta = X/(1-1./alpha);
    std::cout << " Initially, we take prepetual American boundary: theta^0 = " << theta << "\n";
    double A = -1./alpha*pow(theta,1-alpha);
    std::cout << " and: P^0(S) = " << A << " * S^{"<< alpha << "} \n\n";
    
    std::vector<double> yDown(N+1),VDown(N+1);
    { 
      double yMin = log(X*exp(-7.5*sigma*sqrt(dT)));
      double yMax = log(theta);
      double dy = (yMax - yMin)/N;
      // assign grid and corresponding payoff
      for(int i=0;i<=N;i++)
      {
	yDown[i] = yMin + i*dy;
	VDown[i] = X - exp(yDown[i]);
      }
    }
    
    // set y grid up at the final step, using theta^0 as lower limit
    std::vector<double> yUp(N+1),VUp(N+1);
    {
      double yMin = log(theta);
      double yMax = log(X*exp(7.5*sigma*sqrt(dT)));
      double dy = (yMax - yMin)/N;
      // assign grid and corresponding payoff
      for(int i=0;i<=N;i++)
      {
	yUp[i] = yMin + i*dy;        
	VUp[i] = A*exp(alpha*yUp[i]);
      }
    }
    
    std::cout << " Next calculate P^1(S) = A(x)*\int^\\theta B(x,y)*(X-e^y) dy + A(x)*\int_\\theta B(x,y)*V(y,T) dy \n\n";
    // return the value V(x=log(S),0) = A(x)*\int^\theta B(x,y)*(X-e^y) dy + A(x)*\int_\theta B(x,y)*V(y,T) dy
    double value = valueOption(x0,r,sigma,dT,N,yDown,VDown) + valueOption(x0,r,sigma,dT,N,yUp,VUp);
    // get equivalent perpetual value
    if(S<theta)
      std::cout << " P^1(S="<<S<<") = " << value << "  ::  P^0(S="<<S<<") = " << X - S << "\n";
    else
      std::cout << " P^1(S="<<S<<") = " << value << "  ::  P^0(S="<<S<<") = " << A*pow(S,alpha) << "\n";
    
    std::cout << " To find the value of theta^1, we require \n X - theta^1 - P^1(theta^1) = 0 \n In this case we have:- \n";
    std::cout << X << " -  " << S << " - P^1(S="<<S<<") = " << X - S - value << "\n";
  }
  
}