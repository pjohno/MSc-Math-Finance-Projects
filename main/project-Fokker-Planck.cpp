#include "math60082_tridag.hpp"
#include "math60082_lagrange_interp.hpp"
#include <functional>
#include <vector>
#include <cmath>
#include <algorithm>
#include <boost/math/distributions.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
using namespace std;

int main()
{
    double kappa=1,theta=0.,T=10,sigma=0.1;
    double x0=0.,xT=0.;
    
    boost::math::normal OU(x0*exp(-kappa*T) + theta*(1-exp(-kappa*T)),sigma*sqrt((1-exp(-2.*kappa*T))/(2.*kappa)));
    std::cout.precision(8);
    std::cout << " Solution pdf(x=" << xT << ",T="<<T<<") = " << pdf(OU,xT) << endl;
    
    double valueOld=1.;
    for(int k=1;k<=10;k++)
    {
        int n=4*pow(2,k);
        auto start = std::chrono::steady_clock::now(); 
        
        std::vector<double> X(n+1);
        std::vector<double> U(n+1);
        
        double sd_est = sigma*sqrt((1-exp(-2.*kappa*T))/(2.*kappa));
        double xMin=theta-10*sd_est;
        double xMax=theta+10*sd_est;

        int jMax=n;
        int iMax=n;
        double dx=(xMax-xMin)/jMax;
        double dt=T/iMax;
        int jStar0 = (x0-xMin)/dx;
        for(int j=0;j<=n;j++)
        {
            X[j] = xMin+j*dx;
            if(j==jStar0)
                U[j] = 1./dx;
            else
                U[j] = 0.;
        }
        
        for(int i=iMax-1;i>=0;i--)
        {
            // declare vectors for matrix equations
            vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
            // set up matrix equations a[j]=
            a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
            for(int j=1;j<=jMax-1;j++)
            {
                a[j]=0.25*(sigma*sigma/dx/dx + kappa*(theta-X[j])/dx);
                b[j]=-0.5*sigma*sigma/dx/dx + 0.5*kappa - 1./dt;
                c[j]=0.25*(sigma*sigma/dx/dx - kappa*(theta-X[j])/dx);
                d[j]=-a[j]*U[j-1]-(b[j]+2./dt)*U[j]-c[j]*U[j+1];
            }
            a[jMax]= 0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
            MATH60082::thomasSolve(a,b,c,d);
            // set old=new 
            U=d;
        }// finish looping through time steps 
        
        double value=MATH60082::lagrangeInterpolation(U,X,xT);
        double valueExtrap=(4.*value - valueOld)/3.;
        
        auto finish = std::chrono::steady_clock::now(); 
        auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
        
        std::cout << n << " " << value << " " ;
        std::cout << valueExtrap << " ";
        std::cout << " :: ("<< elapsed.count()<< ")"<<std::endl;
        valueOld=value;
                
    }
    
    return 0;
}
