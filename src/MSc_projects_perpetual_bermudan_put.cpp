#include "MSc_projects_perpetual_bermudan_put.hpp"
#include "MSc_projects_quad.hpp"
#include <MSc_projects_newton.hpp>


namespace MSC_PROJECTS
{
   
    double normalDistribution_builtin(double x)
    {
        return 0.5*erfc(-x/sqrt(2.));
    }
    
    double valueBermudanPut(double X,double theta,double  S,double r,double sigma,double dT,int n,
                         const std::vector<double> &y,const std::vector<double> &v)
    {
        double x0 = log(S);
        double alpha = -2*r/sigma/sigma;
        double A = v[n]*exp(-alpha*y[n]);
        double d2=(log(theta/S) - (r-sigma*sigma/2.)*(dT))/(sigma*sqrt(dT));
        double d1=d2 - sigma*sqrt(dT); 
        double d4=(y[n] - x0 - (r-sigma*sigma/2.)*(dT))/(sigma*sqrt(dT));
        double d3=d4 - alpha*sigma*sqrt(dT);
        return X*exp(-r*dT)*normalDistribution_builtin(d2) - S*normalDistribution_builtin(d1) + QUAD::valueOption(x0,r,sigma,dT,n,y,v) + A*pow(S,alpha)*normalDistribution_builtin(-d3);
    }
    
    int valueBermudanPutOption(const std::vector<double> &Si,std::vector<double> &Vi,double &theta,double X,double r,double sigma,double dT,int N,double xi,double tol)
    {
        std::vector<double> y(N+1),vOld(N+1);
        
        // for a perpetual American put option we have 
        double alpha = -2*r/sigma/sigma;
        theta = X/(1-1./alpha);
        double A = -1./alpha*pow(theta,1-alpha);
        
        //theta = X;
        // set y grid up at the final step, using theta^0 as lower limit
        {
            double yMin = log(theta);
            double yMax = log(theta*exp(xi*sigma*sqrt(dT)));
            double dy = (yMax - yMin)/N;
            // assign grid and corresponding payoff
            for(int i=0;i<=N;i++)
            {
                y[i] = yMin + i*dy;        
                vOld[i] = A*exp(alpha*y[i]);
            }
        }
        
        int KMAX=10000;
        for(int k=1;k<=KMAX;k++)
        {
            std::vector<double> x(N+1),vNew(N+1);
            
            // first find theta, using old theta in the integration
            double thetaOld = theta;
            findroot(theta+X/100.,theta,100,tol*X,[&](double S){return X - S - valueBermudanPut(X,thetaOld,S,r,sigma,dT,N,y,vOld);}); 
            
            //         std::cout << k*dT << " S_f = " << theta << std::endl;
            
            { 
                double xMin = log(theta);
                double xMax = log(theta*exp(xi*sigma*sqrt(dT)));
                double dx = (xMax - xMin)/N;
                // assign grid and corresponding payoff
                for(int j=0;j<=N;j++)
                {
                    x[j] = xMin + j*dx;
                    vNew[j] = valueBermudanPut(X,thetaOld,exp(x[j]),r,sigma,dT,N,y,vOld);
                }
            }
            
            // reset 
            y = x;
            vOld = vNew;
            
            if(fabs(thetaOld - theta)<tol*X)break;
            
            if(k == KMAX)throw;
            
        }
        
        std::vector<double> x(Si.size()),vNew(Si.size());
        // assign grid and corresponding payoff
        for(unsigned int j=0;j<Si.size();j++)
        {
            x[j] = log(Si[j]);
            vNew[j] = valueBermudanPut(X,theta,exp(x[j]),r,sigma,dT,N,y,vOld);
        }
        Vi = vNew;
        return 0;
    }
    
}
