#include "MSc_projects_pde_solvers.hpp"
#include <cmath>
using namespace std;

namespace MSC_PROJECTS {
    
    
    // A generic lagrange interpolation function
    double lagrangeInterpolation(const vector<double>& y,const vector<double>& x,double x0,unsigned int n)
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
    )
    {
        int iMax=T*n;
        int jMax=n;
        //set up local parameters
        double dt = (T-t)/iMax;
        double dX = (xMax-xMin)/jMax;
        X.resize(jMax+1);pt.resize(jMax+1);Jt.resize(jMax+1);
        // setup and initialise the stock price 
        for(int j=0;j<=jMax;j++)
        {
            X[j] = xMin + j*dX;
        }
        // setup and initialise the final conditions on the option price 
        // this just uses U(X) as the final condition
        for(int j=0;j<=jMax;j++)
        {
            Jt[j] = pow(X[j],gamma)/gamma;
        }
        
        // start looping through time levels
        for(int i=iMax-1;i>=0;i--)
        {
            vector<double> JtOld(Jt);
            // first calculate optimal p using the values at t+dt
            for(int j=1;j<=jMax-1;j++)
            {
                double dJdX = ( JtOld[j+1] - JtOld[j-1] ) / 2. / dX;
                double d2JdX2 = ( JtOld[j+1] - 2.*JtOld[j] + JtOld[j-1] ) / dX / dX;
		// bound p in between 0 and 5
                double p = std::min(std::max(0.,(r-mu)*X[j]*dJdX / (sigma*sigma*X[j]*X[j]*d2JdX2)),5.);
                if(d2JdX2>0)p=5.;
                pt[j]=p;
            }
            // declare vectors for matrix equations
            vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
            // boundary condition assumes that p=0 at X=0
            a[0]=0.;b[0]=1./dt + .5*r*xMin/dX;c[0]= -0.5*r*xMin/dX;d[0] = JtOld[0]/dt + 0.5*r*xMin*( JtOld[1]-JtOld[0] )/dX;
            for(int j=1;j<=jMax-1;j++)
            {
                double p=pt[j];
                // change these ptOld[j] to p's to update at each step
                if((r+(mu-r)*p)*X[j]/dX>p*p*sigma*sigma*X[j]*X[j]/dX/dX)
                {
                    a[j]=0.5*(p*p*sigma*sigma*X[j]*X[j]/dX/dX);
                    b[j]=-p*p*sigma*sigma*X[j]*X[j]/dX/dX-(r+(mu-r)*p)*X[j]/dX-1./dt;
                    c[j]=0.5*p*p*sigma*sigma*X[j]*X[j]/dX/dX+(r+(mu-r)*p)*X[j]/dX;
                }
                else
                {
                    a[j]=0.5*(p*p*sigma*sigma*X[j]*X[j]/dX/dX-(r+(mu-r)*p)*X[j]/dX);
                    b[j]=-p*p*sigma*sigma*X[j]*X[j]/dX/dX-1./dt;
                    c[j]=0.5*(p*p*sigma*sigma*X[j]*X[j]/dX/dX+(r+(mu-r)*p)*X[j]/dX);
                }
                d[j]=- 1./dt*JtOld[j];
            }
            a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;
            {
                double pHat = (mu - r)/sigma/sigma/(1-gamma);
                double aHat=(r+pHat*(mu-r)) - 0.5*pHat*pHat*sigma*sigma*(1-gamma);
                d[jMax] = pow(X[jMax]*exp(aHat*(T-(i+0.5)*dt)),gamma)/gamma;
            }
            // solve the matrix equations, Jt now stores solution at t+dthalf
            Jt = tridag(a,b,c,d);
        }
        
        return 0;
        
    }// finished time looping
    
    
}
