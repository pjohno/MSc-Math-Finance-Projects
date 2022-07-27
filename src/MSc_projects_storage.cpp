#include "MSc_projects_storage.hpp"
#include "math60082_tridag.hpp"
#include <iostream>
using namespace std;

int MSC_PROJECTS::solveStorageOption(int M,int N,double Qmax,double Xmax,std::vector<double> &Q,std::vector<double> &X,
                                     double r,double p,double sigma,std::vector<std::vector<double>> &v,int qMax,double tol)
{
    Q.resize(M+1);X.resize(2*N+1);
    double dX = Xmax/N;
    double dQ = Qmax/M;
    for(int i=0;i<=2*N;i++)
        X[i] = -Xmax + i*dX;
    for(int j=0;j<=M;j++)
        Q[j] = j*dQ;
    
    v.clear();
    v.resize(2*N+1,std::vector<double>(M+1,0.));
    
    
    
    for(int q=0;q<qMax;q++)
    {
        
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
            MATH60082::thomasSolve(a,b,c,d);
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
            MATH60082::thomasSolve(a,b,c,d);
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
            MATH60082::thomasSolve(a,b,c,d);
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
            MATH60082::thomasSolve(a,b,c,d);
            // update v with solution contained in d
            for(int i=0;i<=N;i++)
                v[N+i][j]=d[i];
        }
        
        double error=0.;
        // to check for errors after sweeping through the regions, only the value at i=N-1 will show up any errors
        for(int j=1;j<=M;j++)
        {
            int i=N-1;
            double kappa_i=sigma*sigma*dQ/dX/dX/2./X[i];
            double rho_i=r*dQ/X[i];
            // calculate residual
            double rj = kappa_i*v[i-1][j] + (1.-2.*kappa_i-rho_i)*v[i][j] + kappa_i*v[i+1][j]  - v[i][j-1]-p*dQ;
            error += rj*rj;
        }
        if(error<tol*tol)break;
    }
    
    return 0;
}
