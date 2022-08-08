#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace std;
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


int main()
{
    // problem parameters
    double gamma_M=-3;
    double gamma_I=-3;
    double a_M=0.1;
    double b_M=0.02;
    double c_M=0.2;
    double r=0.0578;
    double mu=r+0.02;
    
    double mertonRatio=2.;
    
    double sigma=sqrt( (mu - r) / (mertonRatio*(1-gamma_M) ) );
    cout << " Sigma = " << sigma << endl;
    double T=1;
    double H=1;
    double L=0.018;
    double tol=1.e-5;
    double xMin=0.5; // minimum fund value
    double xMax=5; // minimum fund value 
    int n=1000;
    
    int iMax=T*n;
    int jMax=n;
    vector<double> X(jMax+1),pt(jMax+1),Jt(jMax+1),It(jMax+1);
    //set up local parameters
    double dt = T/iMax;
    double dX = (xMax-xMin)/jMax;
    // setup and initialise the stock price 
    for(int j=0;j<=jMax;j++)
    {
        X[j] = xMin + j*dX;
        double Fx=a_M*X[j];
        Jt[j] = pow(Fx,gamma_M)/gamma_M;
        It[j] = pow(X[j]-Fx,gamma_I)/gamma_I;
    }
    
    // start looping through time levels
    for(int i=iMax-1;i>=0;i--)
    {
        vector<double> JtOld(Jt);
        vector<double> ItOld(It);
        
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
        a[0]=0.;b[0]=-1./dt - r*xMin/dX;c[0]= r*xMin/dX;
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
        }
        {
            double pHat = pt[jMax-1];
            double aHat=(r+pHat*(mu-r)) - 0.5*pHat*pHat*sigma*sigma*(1-gamma_M);
            a[jMax]=0.;b[jMax]=-1./dt+gamma_M*aHat;c[jMax]=0.;
        }
        
        for(int j=0;j<=jMax;j++)
            d[j] =-JtOld[j]/dt ;
        // solve the matrix equations, Jt now stores solution at t+dthalf
        thomasSolve(a,b,c,d);
        Jt = d;
        
        // a,b,c stay the same, 
        // except at the boundary
        {
            double pHat = pt[jMax-1];
            double aHat=(r+pHat*(mu-r)) - 0.5*pHat*pHat*sigma*sigma*(1-gamma_I);
            a[jMax]=0.;b[jMax]=-1./dt+gamma_I*aHat;c[jMax]=0.;
        }
        for(int j=0;j<=jMax;j++)
            d[j] = -ItOld[j]/dt;

        thomasSolve(a,b,c,d);
        It = d;
        
    }
    
    std::vector<double> investmentYield_Manager(jMax+1),investmentYield_Investor(jMax+1);
    for(int j=0;j<=jMax;j++)
    {
        // investment yield = 1/T * log(certainty equivalent / initial investment)
        double initialInvest_M=a_M*X[j];
        investmentYield_Manager[j] = 1./T * log(pow( gamma_M * Jt[j] , 1./gamma_M ) / initialInvest_M);
        double initialInvest_I=(1-a_M)*X[j];
        investmentYield_Investor[j] = 1./T * log(pow( gamma_I * It[j] , 1./gamma_I ) / initialInvest_I);
    }
    
    // and output results to file
    ofstream output("test.csv");
    for(int j = 0 ; j<=jMax ; j++)
    {
        string sep = " , ";
        output << X[j] <<sep<< pt[j] <<sep<< Jt[j] <<sep<< It[j] <<sep<< investmentYield_Manager[j] <<sep<< investmentYield_Investor[j] << endl;
    }
    
}

