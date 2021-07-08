#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

// A generic lagrange interpolation function
double lagrangeInterpolation(const vector<double>& y,const vector<double>& x,double x0,unsigned int n=4)
{
    if(x.size()<n)return lagrangeInterpolation(y,x,x0,x.size());
    if(n==0)throw;
    int nHalf = n/2;
    int jStar;
    double dx=x[1]-x[0];
    if(n%2==0)
        jStar = int(x0/dx) -(nHalf-1);
    else
        jStar = int(x0/dx+0.5)-(nHalf);
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

// a tridiagonal solver
vector<double> tridiagonalSolver(const vector<double>& a,const vector<double>& beta,const vector<double>& c,const vector<double>& d)
{
    int n=a.size();
    vector<double> b(n);
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

// payoff for the option, could change this to call option here
double payoff(double S,double X)
{
    return max(X-S,0.);
}

// calculate the option value at t=0 
// output the decision boundary at which we decide to exercise the option
double optionValue(double S0,double X,double T,int nE,double r,double D0,double sigma
,int n,int gridFactor)
{
    ofstream output("test.dat");
    
    int iMax = 2*n*T,jMax = n*gridFactor;
    // nE has to be > 0
    if(nE<=0)throw;
    // declare and initialise local variables (ds,dt)
    double S_max=gridFactor*X;
    double dS=S_max/jMax;
    // now dt is to be the timestep in between exercise dates
    // so divide total number of time steps between exercise periods
    int iSteps = max(1,iMax/nE);
    // dt gives variation in between exercise dates
    double dt=(T/nE)/iSteps;
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
        double alpha = -2*r/sigma/sigma;
        double theta = X/(1-1./alpha);
        double A = -1./alpha*pow(theta,1-alpha);
        double payoff ;
        if(S[j] < theta) payoff = X-S[j];
        else payoff = A*pow(S[j],alpha);
        vOld[j] = payoff;
        vNew[j] = payoff;
    }
    
    // start looping through exercise dates
    for(int exDate=nE-1;exDate>=0;exDate--)
    {
        vector<int> exerciseDecision(jMax+1); 
        for(int i=iSteps-1;i>=0;i--)
        {
            // declare vectors for matrix equations
            vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
            // set up matrix equations a[j]=
            a[0]=0.;b[0]=1.;c[0]=0.;d[0] = X*exp(-r*(T-i*dt));
            for(int j=1;j<=jMax-1;j++)
            {
                a[j]=0.25*(sigma*sigma*j*j-(r-D0)*j);
                b[j]=-0.5*sigma*sigma*j*j - 0.5*r - 1./dt;
                c[j]=0.25*(sigma*sigma*j*j+(r-D0)*j);
                d[j]=-a[j]*vOld[j-1]-(b[j]+2./dt)*vOld[j]-c[j]*vOld[j+1];
            }
            double alpha = -2*r/sigma/sigma;
            a[jMax]= -jMax - 0.5*alpha;b[jMax]=jMax - 0.5*alpha;c[jMax]=0.;d[jMax] = 0.;
            vNew = tridiagonalSolver(a,b,c,d);
            // set old=new 
            vOld=vNew;
        }// finish looping through time steps for a single exercise period
        
        
        // now apply exercise conditions
        for(int j=0;j<=jMax;j++)
        {
            if(vNew[j]<payoff(S[j],X))
            {
                vNew[j] = payoff(S[j],X);
                exerciseDecision[j]=1;
            }
            else
                exerciseDecision[j]=0;
        }
        // record exercise decision 
        int current=exerciseDecision[0];
        for(int j=0;j<=jMax;j++)
        {
            if(exerciseDecision[j]!=current)
            {
                double boundaryprice = (-vOld[j-1] + (vOld[j] - vOld[j - 1]) / (S[j] - S[j - 1])*S[j - 1] + X) / (1 + (vOld[j] - vOld[j - 1]) / (S[j] - S[j - 1]));
                output << " " << T- exDate * (T/nE) << " , " << boundaryprice << endl; 
                cout << " Estimated exercise boundary at time t= " << exDate * (T/nE) << " is " << boundaryprice << endl; 
                break;
            }
            current = exerciseDecision[j];
        }
        // set old=new 
        vOld=vNew;
    } // go through all exercise dates
    return lagrangeInterpolation(vNew,S,S0);
    
}

/* 
 * Template code for solving a Bermudan put option
 */
int main()
{
    double S0=100.,X=100.,r=0.1,d=0.,sigma=0.5;
    
    // and the length of the contract
    double T=100,DT=0.1;
    // try varying number of exercise dates
    int exerciseDates=T/DT;
    
    // numerical grid parameters, increase n to get more accuracy
    int n=400,gridFactor=5;
    
    cout << " Running code with n=" << n << endl;
    double value=optionValue(S0,X,T,exerciseDates,r,d,sigma,n,gridFactor);
    cout << " Value of the a Bermudan put option with "<< exerciseDates << " exercise dates ";
    cout << "V(S="<<S0<<",t=0;T="<<T <<") := " << value << endl;
    
}
/*#######################
 *       OUTPUT
 *#######################
 *  Running code with n=100
 * Estimated exercise boundary at time t=0.8 is S_f=0.865
 * Estimated exercise boundary at time t=0.6 is S_f=0.815
 * Estimated exercise boundary at time t=0.4 is S_f=0.785
 * Estimated exercise boundary at time t=0.2 is S_f=0.765
 * Estimated exercise boundary at time t=0 is S_f=0.745
 * Value of the a Bermudan put option with 5 exercise dates V(S=1,t=0;T=1) := 0.0975024
 * Running code with n=400
 * Estimated exercise boundary at time t=0.8 is S_f=0.86375
 * Estimated exercise boundary at time t=0.6 is S_f=0.81625
 * Estimated exercise boundary at time t=0.4 is S_f=0.78375
 * Estimated exercise boundary at time t=0.2 is S_f=0.76375
 * Estimated exercise boundary at time t=0 is S_f=0.74625
 * Value of the a Bermudan put option with 5 exercise dates V(S=1,t=0;T=1) := 0.0975141
 * Running code with n=1600
 * Estimated exercise boundary at time t=0.8 is S_f=0.863438
 * Estimated exercise boundary at time t=0.6 is S_f=0.815313
 * Estimated exercise boundary at time t=0.4 is S_f=0.784687
 * Estimated exercise boundary at time t=0.2 is S_f=0.762813
 * Estimated exercise boundary at time t=0 is S_f=0.745313
 * Value of the a Bermudan put option with 5 exercise dates V(S=1,t=0;T=1) := 0.097515
 * 
 */
