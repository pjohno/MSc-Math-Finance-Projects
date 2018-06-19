#include "MSc_projects_quad.hpp"
#include "MSc_projects_table.hpp"
#include <iostream>
#include <iomanip>
using namespace MSC_PROJECTS;

double callOption(double  S,double X,double r,double sigma,double T,int N=500,double xi=7.5)
{
    // set y grid up at the final step, using strike price as lower limit
    std::vector<double> y(N+1),V(N+1);
    QUAD::resetGrid(N,log(X),log(X*exp(xi*sigma*sqrt(T))),y);
    // assign payoff, since e^y>=X we don't need max function
    for(unsigned int i=0;i<y.size();i++)
        V[i] = exp(y[i])-X;
    
    // return the value V(x=log(S),0) = A(x)*\int B(x,y)*V(y,T) dy
    return QUAD::valueOption(S,r,sigma,T,N,y,V);
}

double upAndOutBarrierCallOption(double  S,double X,double r,double sigma,double T,double B,int N)
{
    double dT = T/2.;
    
    // storage for the solution
    std::vector<double> x(N+1),y(N+1),vOld(N+1),vNew(N+1);
    
    // set x grid up with the barrier at the top
    QUAD::resetGrid(N,log(S*exp(-7.5*sigma*sqrt(dT))), log(B),x);
    // set y grid up at the final step, using strike price as lower limit
    QUAD::resetGrid(N,log(X),log(X*exp(7.5*sigma*sqrt(T))),y);
    
    // assign payoff, since e^y>=X we don't need max function
    for(unsigned int i=0;i<y.size();i++)
        vOld[i] = exp(y[i])-X;
    
    // get the value V(x_i,T/2) = A(x)*\int B(x,y)*V(y,T) dy
    for(unsigned int i=0;i<x.size();i++)
        vNew[i] = QUAD::valueOption(exp(x[i]),r,sigma,dT,N,y,vOld);
    
    // overwrite values of y and vOld to be 
    y = x;
    vOld = vNew;
    // return the value V(x=log(S),0) = A(x)*\int B(x,y)*V(y,T/2) dy
    return QUAD::valueOption(S,r,sigma,dT,N,y,vOld);
    
}

int main()
{
    double  S = 100., X = 100., r = 0.06, sigma = 0.2, T = 1.;
    std::cout << " A call option using QUAD...\n\n";
    tableRow("N","V(S,0;N)","R","conv rate");
    emptyTableRow(4);
    double valueOld = 1.,diffOld=1.;
    int k=2;
    for(int n=10;n<=5000;n*=k)
    {
        double value = callOption(S,X,r,sigma,T,n);
        double diff = value - valueOld;
        double R = diffOld/diff;
        tableRow( n ,value ,R ,log(R)/log(k));
        valueOld = value;
        diffOld = diff;
    }
    
    std::cout << "\n\n A barrier call option using QUAD...\n";
    tableRow("N","V(S,0;N)","R","conv rate");
    emptyTableRow(4);
    for(int n=10;n<=5000;n*=2)
    {
        double value = upAndOutBarrierCallOption(S,X,r,sigma,T,110.,n);
        double diff = value - valueOld;
        double R = diffOld/diff;
        tableRow( n ,value ,R ,log(R)/log(k));
        valueOld = value;
        diffOld = diff;
    }
}
