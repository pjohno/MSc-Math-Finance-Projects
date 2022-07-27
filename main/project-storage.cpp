#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>
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
    /* Template code for the storage valuation
     *  
     *      X dV/dQ + 1/2 sigma^2 d^2V/dX^2 -rV - p min(X,0) = 0
     *  
     *  At X->-\infty we have
     *     V = pQ
     *  As X->infty we have
     *     V = 0.
     * 
     */
    // load parameters and storage
    int M=10;
    int N=25;
    double Qmax=1.;
    double Xmax=25.;
    std::vector<double> Q(M+1),X(2*N+1);
    
    double dX = Xmax/N;
    double dQ = Qmax/M;
    for(int i=0;i<=2*N;i++)
        X[i] = -Xmax + i*dX;
    for(int j=0;j<=M;j++)
        Q[j] = j*dQ;
    
    double r=0.01;
    double p=1;
    double sigma=0.5;
    
    std::vector<std::vector<double>> v(2*N+1,std::vector<double>(M+1,0.));
    // SOLVE FOR REGION A
    
    // solve at j=0 according to (4) and (5)
    {
        // store tridiagonal system
        vector<double> a(N+1),b(N+1),c(N+1),d(N+1);
        // constant kappa after rearranging (4)
        double kappa=2.*r*dX*dX/sigma/sigma;
        
        a[0]=0.;b[0]=1.;c[0]=0.;d[0]=0.;
        for(int i=1;i<N;i++)
        {
            a[i]=1.;b[i]=-2.-kappa;c[i]=1;d[i]=0;
        }
        a[N]=1.;b[N]=-2.-kappa;c[N]=0.;d[N]=-v[N+1][0];
        
        // solve
        thomasSolve(a,b,c,d);
        // update v with solution
        for(int i=0;i<=N;i++)
            v[i][0]=d[i];
    }
    // solve at j=1 to M according to (7) and (8)
    for(int j=1;j<=M;j++)
    {
        // store tridiagonal system
        vector<double> a(N+1),b(N+1),c(N+1),d(N+1);
        
        a[0]=0.;b[0]=1.;c[0]=0.;d[0]=Q[j];
        for(int i=1;i<N;i++)
        {
            double kappa_i=sigma*sigma*dQ/dX/dX/2./X[i];
            double rho_i=r*dQ/X[i];
            a[i]=kappa_i;b[i]=1.-2.*kappa_i-rho_i;c[i]=kappa_i;d[i]=v[i][j-1]+p*dQ;
        }
        // constant kappa after rearranging 
        double kappa=2.*r*dX*dX/sigma/sigma;
        a[N]=1.;b[N]=-2.-kappa;c[N]=0.;d[N]=-v[N+1][j];
        
        // solve
        thomasSolve(a,b,c,d);
        // update v with solution contained in d
        for(int i=0;i<=N;i++)
            v[i][j]=d[i];
    }
    
    // SOLVE FOR REGION B
    
    // solve at j=M according to (9) and (10)
    {
        // store tridiagonal system
        vector<double> a(N+1),b(N+1),c(N+1),d(N+1);
        // constant kappa after rearranging (4)
        double kappa=2.*r*dX*dX/sigma/sigma;
        
        a[0]=0.;b[0]=-2.-kappa;c[0]=1.;d[0]=-v[N-1][M];
        for(int i=1;i<N;i++)
        {
            a[i]=1.;b[i]=-2.-kappa;c[i]=1;d[i]=0;
        }
        a[N]=0.;b[N]=1.;c[N]=0.;d[N]=0.;
        
        // solve
        thomasSolve(a,b,c,d);
        // update v with solution
        for(int i=0;i<=N;i++)
            v[N+i][M]=d[i];
    }
    // solve at j=M-1 to 0 according to (12) and (13)
    for(int j=M-1;j>=0;j--)
    {
        // store tridiagonal system
        vector<double> a(N+1),b(N+1),c(N+1),d(N+1);
        
        // constant kappa after rearranging 
        double kappa=2.*r*dX*dX/sigma/sigma;
        a[0]=0.;b[0]=-2.-kappa;c[0]=1.;d[0]=-v[N-1][j];
        for(int i=1;i<N;i++)
        {
            double kappa_i=sigma*sigma*dQ/dX/dX/2./X[N+i];
            double rho_i=r*dQ/X[N+i];
            a[i]=kappa_i;b[i]=-1.-2.*kappa_i-rho_i;c[i]=kappa_i;d[i]=-v[N+i][j+1];
        }
        a[N]=0.;b[N]=1.;c[N]=0.;d[N]=0.;
        
        // solve
        thomasSolve(a,b,c,d);
        // update v with solution contained in d
        for(int i=0;i<=N;i++)
            v[N+i][j]=d[i];
    }
    
    
    
}

