#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace std;

vector<double> tridag(const vector<double>& a,const vector<double>& beta,const vector<double>& c,vector<double>& d)
{
    int n=a.size();
    vector<double> b(n);
    // move d to rhs
    vector<double> rhs(d);
    // initial first value of b
    b[0]=beta[0];
    for(int j=1;j<n;j++)
    {
        b[j]=beta[j]-c[j-1]*a[j]/b[j-1];  
        rhs[j]=rhs[j]-rhs[j-1]*a[j]/b[j-1];
    }
    // calculate solution
    rhs[n-1]=rhs[n-1]/b[n-1];
    for(int j=n-2;j>=0;j--)
        rhs[j]=(rhs[j]-c[j]*rhs[j+1])/b[j];
    return rhs;
}

// function to calculate jump integral term
double trapeziumRule(double Sj,const std::vector<double> &S,
                     const std::vector<double> &v,
                     double J_gamma,double J_mu)
{
    double sum=0.;
    for(int k=0;k<S.size()-1;k++)
    {
        double dS = (S[k+1]-S[k]);
        // value of S at mid point
        double Smid = 0.5*(S[k+1]+S[k]);
        double mu = log(Smid/Sj)-J_mu;
        double sig = 2*J_gamma*J_gamma;
        sum += dS*(v[k+1]+v[k])*exp((-mu*mu)/sig)/(sqrt(sig*M_PI)*Smid);
    }
    return 0.5*sum;
}

int main()
{
    // declare and initialise Black Scholes parameters
    double S0=100,X=100.,T=0.25,r=0.05,sigma=0.15;
    // declare jump stuff
    double J_lambda=0.1,J_gamma=0.45;
    double J_mu = -0.9;
    
    // declare and initialise grid paramaters 
    int iMax=100,jMax=100;
    // declare and initialise local variables (ds,dt)
    double S_max=5*X;
    double dS=S_max/jMax;
    double dt=T/iMax;
    // create storage for the stock price and option price (old and new)
    vector<double> S(jMax+1),vOld(jMax+1),vNew(jMax+1);
    // setup and initialise the stock price 
    for(int j=0;j<=jMax;j++)
    {
        S[j] = j*dS;
    }
    // setup and initialise the final conditions on the option price 
    for(int j=0;j<=jMax;j++)
    {
        vOld[j] = max(X-S[j],0.);
        vNew[j] = max(X-S[j],0.);
    }
    // start looping through time levels
    for(int i=iMax-1;i>=0;i--)
    {
        // declare vectors for matrix equations
        vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
        // set up matrix equations a[j]=
        a[0] = 0.;b[0] = 1.;c[0] = 0.;
        d[0] = X*exp(-r*(iMax-i)*dt);
        for(int j=1;j<jMax;j++)
        {
            // calculate jump components
            double jumpDrift = J_lambda*(exp(J_mu+J_gamma*J_gamma/2.)-1.);
            double J_integral_Term = trapeziumRule(S[j],S,vOld,J_gamma,J_mu);
            
            // set up scheme
            a[j] = 0.5*(sigma*sigma*j*j-(r-jumpDrift)*j);
            b[j] =-sigma*sigma*j*j - (r+J_lambda) - 1./dt;
            c[j] = 0.5*(sigma*sigma*j*j+(r-jumpDrift)*j);
            
            // set up right hand side
            d[j] = -vOld[j]/dt - J_lambda * J_integral_Term;
        }
        a[jMax] = 0.;b[jMax] = 1.;c[jMax] = 0.;
        d[jMax] = 0.;
        // solve with Thomas 
        vNew = tridag(a,b,c,d);
        // set old=new
        vOld=vNew;
    }
    // finish looping through time levels
    
    // output the estimated option price
    double optionValue;
    {
        int jStar=S0/dS;
        double sum=0.;
        sum+=(S0 - S[jStar])/dS * vNew[jStar+1];
        sum+=(S[jStar+1] - S0)/dS * vNew[jStar];
        optionValue = sum;
    }
    
    cout << " V:= " <<  optionValue << "\n";
    
}
