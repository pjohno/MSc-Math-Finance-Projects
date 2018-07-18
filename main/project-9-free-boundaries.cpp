#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>
using namespace std;

class PerpetualPutOption
{
  
    double X,A,alpha,Sf;
    
public:
    
    PerpetualPutOption(double X,double r,double D0,double sigma):X(X)
    {
        double rDhalfsig2 = r-D0-0.5*sigma*sigma;
        double d = rDhalfsig2*rDhalfsig2 + 2*r*sigma*sigma;
        alpha = (-rDhalfsig2 - sqrt(d))/sigma/sigma;
        Sf = X/(1-1./alpha);
        A = -1./alpha*pow(Sf,1-alpha);
    }
    
    double operator()(double S) const { 
        if(S < Sf)
            return X - S;
        else
            return A*pow(S,alpha);
    }
    
    double getSf() const {return Sf;}
    
};

class PerpetualCallOption
{
  
    double X,A,alpha,Sf;
    
public:
    
    PerpetualCallOption(double X,double r,double D0,double sigma):X(X)
    {
        double rDhalfsig2 = r-D0-0.5*sigma*sigma;
        double d = rDhalfsig2*rDhalfsig2 + 2*r*sigma*sigma;
        alpha = (-rDhalfsig2 + sqrt(d))/sigma/sigma;
        Sf = X/(1-1./alpha);
        A = 1./alpha*pow(Sf,1-alpha);
    }
    
    double operator()(double S) const { 
        if(S > Sf)
            return S - X;
        else
            return A*pow(S,alpha);
    }
    
    double getSf() const {return Sf;}
    
};

int main()
{
    PerpetualPutOption P(100,0.06,0.0,0.2);
    PerpetualCallOption C(100,0.06,0.0,0.2);
    for(int i = 0;i<=100;i++)
    {
        double S = i*2;
        cout << S << " " << P(S) << " " << C(S) << endl;
    }
}
