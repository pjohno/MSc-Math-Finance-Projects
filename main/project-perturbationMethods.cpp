#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

double normalDistribution(double x)
{
    return 0.5*erfc(-x/sqrt(2.));
}

// return the value of a put option using the black scholes formula
double putOptionPrice(double S,double T,double X,double r,double sigma)
{
    if(S<1.e-14)return X*exp(-r*T); // check if asset worthless
    if(sigma<1.e-14) // check if sigma zero
    {
        if(S<X*exp(-r*T))return X*exp(-r*T)-S;
        else return 0.;
    }
    if(fabs(T)<1.e-14) // check if we are at maturity
    {
        if(S<X)return X-S;
        else return 0.;
    }
    // calculate option price
    double d1=(log(S/X) + (r+sigma*sigma/2.)*T)/(sigma*sqrt(T));
    double d2=(log(S/X) + (r-sigma*sigma/2.)*T)/(sigma*sqrt(T));
    return  normalDistribution(-d2)*X*exp(-r*T) - normalDistribution(-d1)*S;
}
/* 
 *    ON INPUT:
 *    a, b and c -- are the diagonals of the matrix
 *    rhs        -- is the right hand side
 *    ON OUTPUT:
 *    a, b, c             -- unchanged
 *    rhs                 -- solution to Ax=b
 */
void thomasSolve(const std::vector<double> &a,const std::vector<double> &b_,const std::vector<double> &c,std::vector<double> &rhs)
{
    int n=a.size();
    std::vector<double> b(n);
    // initial first value of b
    b[0]=b_[0];
    for(int j=1;j<n;j++)
    {
        b[j]=b_[j]-c[j-1]*a[j]/b[j-1];  
        rhs[j]=rhs[j]-rhs[j-1]*a[j]/b[j-1];
    }
    // calculate solution
    rhs[n-1]=rhs[n-1]/b[n-1];
    for(int j=n-2;j>=0;j--)
        rhs[j]=(rhs[j]-c[j]*rhs[j+1])/b[j];
}
// A generic lagrange interpolation function
double lagrangeInterpolation(const vector<double>& y,const vector<double>& x,double x0,unsigned int n=4)
{
    if(x.size()<n)return lagrangeInterpolation(y,x,x0,x.size());
    if(n==0)throw;
    int nHalf = n/2;
    int jStar;
    double dx=x[1]-x[0];
    if(n%2==0)
        jStar = int((x0 - x[0])/dx) -(nHalf-1);
    else
        jStar = int((x0 - x[0])/dx+0.5)-(nHalf);
    jStar=std::max(0,jStar);
    jStar=std::min(int(x.size()-n),jStar);
    if(n==1)return y[jStar];
    double temp = 0.;
    for(unsigned int i=jStar;i<jStar+n;i++){
        double  int_temp;
        int_temp = y[i];
        for(unsigned int j=jStar;j<jStar+n;j++){
            if(j==i){continue;}
            int_temp *= ( x0 - x[j] )/( x[i] - x[j] );
        }
        temp += int_temp;
    }
    // end of interpolate
    return temp;
}

void caluculateOptionPrice( 
vector<double> &sHat,
vector<vector<double>> &value,
double T, // contract maturity date
double X, // risk free interest rate
double r, // risk free interest rate
// declare number of terms in the approximation
int n,
// declare and initialise grid paramaters 
int iMax, // number of time steps
int jMax, // number of space steps
int upperBoundary=5 // place upper boundary at this integer times the strike price -- default to 5
)
{
    
    // declare and initialise local grid variables (dlambda,dt)
    double sHat_min=-5*X;
    double sHat_max=5*X;
    double dS=(sHat_max-sHat_min)/jMax;
    double dt=T/iMax;
    
    std::vector<std::vector<double>> vNew(n+1,std::vector<double>(jMax+1));
    // setup and initialise the final conditions on the option price 
    // for ni=0 we have P_0 = max(-\hat S,0) 
    for(int j=0;j<=jMax;j++)
    {
        sHat[j] = sHat_min + j*dS;;
    }
    for(int j=0;j<=jMax;j++)
    {
        vNew[0][j] = std::max(-sHat[j],0.);
    }
    // for ni>0 we have P_ni = 0
    for(int ni=1;ni<=n;ni++)
    {
        for(int j=0;j<=jMax;j++)
        {
            vNew[ni][j] = 0.;
        }
    }
    
    for(int i=iMax-1;i>=0;i--)
    {
        std::vector<std::vector<double>> vOld = vNew;
        
        double tau=T-(i+0.5)*dt;
        // declare vectors for matrix equations
        std::vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
        
        // O(sigma) equation
        
        // set up matrix equations for the boundary at lambda=0
        a[0]=0.;b[0]=-1/dS;c[0]=1./dS;d[0] = -1.;
        // set up the scheme in the middle
        for(int j=1;j<=jMax-1;j++)
        {
            a[j]=0.25*X*X*exp(-2.*r*tau)/dS/dS-0.25*r*sHat[j]/dS;
            b[j]=-0.5*X*X*exp(-2.*r*tau)/dS/dS - 0.5*r - 1./dt;
            c[j]=0.25*X*X*exp(-2.*r*tau)/dS/dS+0.25*r*sHat[j]/dS;
            d[j]=-a[j]*vOld[0][j-1]-(b[j]+2./dt)*vOld[0][j]-c[j]*vOld[0][j+1];
        }
        // set up boundary at lambda_max
        a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
        // solve the system of equations with Thomas algorithm
        // note that "d" contains the solution on exit
        thomasSolve(a,b,c,d);
        vNew[0]=d;
        
        
        // O(sigma^2) equation
        
        // set up matrix equations for the boundary at lambda=0
        a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
        // set up the scheme in the middle
        for(int j=1;j<=jMax-1;j++)
        {
            // a,b,c are the same at every order, so don't need to recalculate
            // get the value of d2P_0dS^2 at the grid point, and add it in
            double d2P0dS2 = 0.5*(vNew[0][j-1]-2*vNew[0][j]+vNew[0][j+1] + vOld[0][j-1]-2*vOld[0][j]+vOld[0][j+1])/dS/dS;
            d[j]=-a[j]*vOld[1][j-1]-(b[j]+2./dt)*vOld[1][j]-c[j]*vOld[1][j+1] - sHat[j]*X*exp(-r*tau)*d2P0dS2;
        }
        // set up boundary at lambda_max
        a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
        // solve the system of equations with Thomas algorithm
        // note that "d" contains the solution on exit
        thomasSolve(a,b,c,d);
        vNew[1]=d;
        
        // O(sigma^{n+1}) for n>=2
        // can you fill this in?
        for(int ni=2;ni<=n;ni++)
        {
            // set up matrix equations for the boundary at lambda=0
            a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
            // set up the scheme in the middle
            for(int j=1;j<=jMax-1;j++)
            {
                // a,b,c are the same at every order, so don't need to recalculate
                // what should d be?
                d[j] = 0.;
            }
            // set up boundary at lambda_max
            a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
            // solve the system of equations with Thomas algorithm
            // note that "d" contains the solution on exit
            thomasSolve(a,b,c,d);
            vNew[ni]=d;
            
        }
    }
    
    value = vNew;
    
}

double getOptionValue(const std::vector<double> &sHat,
                      const std::vector<std::vector<double>>& v,
                      double S0, // initial stock value
                      double T, // contract maturity date
                      double X, // risk free interest rate
                      double r, // risk free interest rate
                      double sigma, // volatility of default rate
                      unsigned int n // nth term expansion
)
{
    // finish looping through time levels
    // return the interpolated value at lambda0
    double sHat0 = (S0-X*exp(-r*T))/sigma;
    // 
    if(sHat0 > sHat[sHat.size()-1])
        return 0.;
    else if(sHat0 < sHat.front())
        return X*exp(-r*T)-S0;
    else
    {
        double temp=0.;
        for(int ni=0;ni<=std::min(n,(unsigned int)(v.size())-1);ni++){
            temp+=pow(sigma,ni+1)*lagrangeInterpolation(v[ni],sHat,sHat0);
        }
        return temp;
    }
} 

/* Template code for the Crank Nicolson Finite Difference
 */
int main()
{
    // parameters for the problem
    double T=1,sigma=0.1,r=0.03,X=1.;
    // n is the maximum number of terms to calculate
    // iMax is the number of time steps
    // jMax is the number of space steps
    // sMax is related to the domain of the grid
    int n=1,iMax=500,jMax=500,sMax=5;
    
    vector<double> sHat(jMax+1);
    vector<vector<double>> value(n+1,vector<double>(jMax+1));
    
    // calculate all the solutions P_i for 0<=i<=n
    caluculateOptionPrice(sHat,value,T,X,r,n,iMax,jMax,sMax);
    
    vector<double> errors(n+1,0.);
    
    // print out solutions at different values of S
    ofstream output("test.csv");
    output.precision(12);
    for(int i=-100;i<=100;i++)
    {
        double S0 = X*(1+i/100.);
        double analyticValue = putOptionPrice(S0,T,X,r,sigma);
        output << S0;
        for(int ni=0;ni<=n;ni++)
        {
            double P_ni_Solution = getOptionValue(sHat,value,S0,T,X,r,sigma,ni);
            output << " , " << P_ni_Solution;
            errors[ni] = max(errors[ni],fabs(P_ni_Solution-analyticValue));
        }
        output << " , " << analyticValue << endl;
    }
    
    // create a table with the truncation errors
    cout << " Output written to file \"test.csv\". Truncation errors of the series as follows::\n C_i \t| Error \n";
    for(int ni=0;ni<=n;ni++)
    {
        cout << ni << " \t| " << errors[ni] << "\n";
    }
    cout <<endl;
}
