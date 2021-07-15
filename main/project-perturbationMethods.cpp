#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace std;

double normalDistribution(double x)
{
  static const double RT2PI = sqrt(4.0*acos(0.0));
  static const double SPLIT = 10./sqrt(2);
  static const double a[] = {220.206867912376,221.213596169931,112.079291497871,33.912866078383,6.37396220353165,0.700383064443688,3.52624965998911e-02};
  static const double b[] = {440.413735824752,793.826512519948,637.333633378831,296.564248779674,86.7807322029461,16.064177579207,1.75566716318264,8.83883476483184e-02};
  
  const double z = fabs(x);
  double Nz = 0.0;
  
  // if z outside these limits then value effectively 0 or 1 for machine precision
  if(z<=37.0)
  {
    // NDash = N'(z) * sqrt{2\pi}
    const double NDash = exp(-z*z/2.0)/RT2PI;
    if(z<SPLIT)
    {
      const double Pz = (((((a[6]*z + a[5])*z + a[4])*z + a[3])*z + a[2])*z + a[1])*z + a[0];
      const double Qz = ((((((b[7]*z + b[6])*z + b[5])*z + b[4])*z + b[3])*z + b[2])*z + b[1])*z + b[0];
      Nz = RT2PI*NDash*Pz/Qz;
    }
    else
    {
      const double F4z = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
      Nz = NDash/F4z;
    }
  }
  return x>=0.0 ? 1-Nz : Nz;
}

// return the value of a put option using the black scholes formula
double putOptionPrice(double S,double t,double X,double r,double d,double sigma,double T)
{
  if(S<1.e-14)return X*exp(-r*(T-t)); // check if asset worthless
  if(sigma<1.e-14) // check if sigma zero
  {
    if(S*exp(-d*(T-t))<X*exp(-r*(T-t)))return X*exp(-r*(T-t))-S*exp(-d*(T-t));
    else return 0.;
  }
  if(fabs(T-t)<1.e-14) // check if we are at maturity
  {
    if(S<X)return X-S;
    else return 0.;
  }
  // calculate option price
  double d1=(log(S/X) + (r-d+sigma*sigma/2.)*(T-t))/(sigma*sqrt(T-t));
  double d2=(log(S/X) + (r-d-sigma*sigma/2.)*(T-t))/(sigma*sqrt(T-t));
  return  normalDistribution(-d2)*X*exp(-r*(T-t)) - normalDistribution(-d1)*S*exp(-d*(T-t));
}

// a tridiagonal solver, where a, b and c are the diagonals from the matrix A and d is the right hand side
// return the solution vector
vector<double> tridiagonalSolver(const vector<double>& a,const vector<double>& beta,const vector<double>& c,vector<double>& d)
{
  int n=a.size();
  vector<double> b(n);
  // move d to solution vector
  vector<double> sol(d);
  // initial first value of b
  b[0]=beta[0];
  for(int j=1;j<n;j++)
  {
    b[j]=beta[j]-c[j-1]*a[j]/b[j-1];  
    sol[j]=sol[j]-sol[j-1]*a[j]/b[j-1];
  }
  // calculate solution
  sol[n-1]=sol[n-1]/b[n-1];
  for(int j=n-2;j>=0;j--)
    sol[j]=(sol[j]-c[j]*sol[j+1])/b[j];
  return sol;
}

// for a set of data point (x_i,y_i), use a cubic lagrange interpolation to calculate the value of y(x0)
double interpolate(double x0,double dx,const vector<double> &x,const vector<double> &y)
{
  int jStar= int((x0-x[0])/dx);
  jStar=std::max(1,jStar);
  jStar=std::min(int(x.size()-3),jStar);
  double temp = 0.;
  for(int i=jStar-1;i<jStar+3;i++){
    double  int_temp;
    int_temp = y[i];
    for(int j=jStar-1;j<jStar+3;j++){
      if(j==i){continue;}
      int_temp *= ( x0 - x[j] )/( x[i] - x[j] );
    }
    temp += int_temp;
  }
  // end of interpolate
  return temp;
}

double outerSolution(double S,double t,double T,double X,double r,double d)
{
  if(S*exp(-d*(T-t))>X*exp(-r*(T-t)))
    return X*exp(-r*(T-t))-S*exp(-d*(T-t));
  else
    return 0.;
}

void caluculateOptionPrice( 
vector<double> &sHat,
vector<vector<double>> &value,
double t0, // initial time
double T, // contract maturity date
double r, // risk free interest rate
double X, // risk free interest rate
// declare number of terms in the approximation
int n,
// declare and initialise grid paramaters 
int iMax, // number of time steps
int jMax, // number of space steps
int upperBoundary=5 // place upper boundary at this integer times the strike price -- default to 5
)
{
  // declare and initialise local grid variables (dlambda,dt)
  double sHat_min=-upperBoundary*X;
  double sHat_max=upperBoundary*X;
  double dShat=(sHat_max-sHat_min)/jMax;
  double dt=(T-t0)/iMax;
  // create storage for the default value and derivative price (old and new)
  sHat.resize(jMax+1);
  vector<vector<double>> vOld(n+1,vector<double>(jMax+1)),vNew(n+1,vector<double>(jMax+1));
  // setup and initialise shat
  for(int j=0;j<=jMax;j++)
  {
    sHat[j] = sHat_min + j*dShat;
  }
  // setup and initialise the final conditions on the option price 
  // for n=0 we have P_0 = max(-\hat S,0) 
  for(int j=0;j<=jMax;j++)
  {
    // 
    vNew[0][j] = max(-sHat[j],0.);
  }
  for(int ni=1;ni<=n;ni++)
  {
    for(int j=0;j<=jMax;j++)
    {
      // 
      vNew[ni][j] = 0.;
    }
  }
  
  // start looping through time levels
  for(int i=iMax-1;i>=0;i--)
  {
    double t=t0+(i+0.5)*dt;
    double tau=T-t;
    vOld=vNew;
    // declare vectors for matrix equations
    vector<double> a(jMax+1),b(jMax+1),c(jMax+1),d(jMax+1);
    
    // O(sigma) equation
    
    // set up matrix equations for the boundary at lambda=0
    a[0]=0.;b[0]=-1/dShat;c[0]=1./dShat;d[0] = -1.;
    // set up the scheme in the middle
    for(int j=1;j<=jMax-1;j++)
    {
      a[j]=0.25*X*X*exp(-2.*r*tau)/dShat/dShat-0.25*r*sHat[j]/dShat;
      b[j]=-0.5*X*X*exp(-2.*r*tau)/dShat/dShat - 0.5*r - 1./dt;
      c[j]=0.25*X*X*exp(-2.*r*tau)/dShat/dShat+0.25*r*sHat[j]/dShat;
      d[j]=-a[j]*vOld[0][j-1]-(b[j]+2./dt)*vOld[0][j]-c[j]*vOld[0][j+1];
    }
    // set up boundary at lambda_max
    a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
    // solve the system of equations with Thomas algorithm
    vNew[0] = tridiagonalSolver(a,b,c,d);
    
    // O(sigma^2) equation
    
    // set up matrix equations for the boundary at lambda=0
    a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
    // set up the scheme in the middle
    for(int j=1;j<=jMax-1;j++)
    {
      double d2P0dS2 = 0.5*(vNew[0][j-1]-2*vNew[0][j]+vNew[0][j+1] + vOld[0][j-1]-2*vOld[0][j]+vOld[0][j+1])/dShat/dShat;
      d[j]=-a[j]*vOld[1][j-1]-(b[j]+2./dt)*vOld[1][j]-c[j]*vOld[1][j+1] - sHat[j]*X*exp(-r*tau)*d2P0dS2;
    }
    // set up boundary at lambda_max
    a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
    // solve the system of equations with Thomas algorithm
    vNew[1] = tridiagonalSolver(a,b,c,d);
    
    // O(sigma^n) equation
    
    for(int ni=2;ni<=n;ni++)
    {
      // set up matrix equations for the boundary at lambda=0
      a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
      // set up the scheme in the middle
      for(int j=1;j<=jMax-1;j++)
      {
          /***********
           * 
           * 
           * CAN YOU FILL THIS IN???
           * 
           * 
           */
          a[0]=0.;b[0]=1.;c[0]=0.;d[0] = 0.;
      }
      // set up boundary at lambda_max
      a[jMax]=0.;b[jMax]=1.;c[jMax]=0.;d[jMax] = 0.;
      // solve the system of equations with Thomas algorithm
      // here ni=n-1 
      vNew[ni] = tridiagonalSolver(a,b,c,d);
    }
  }
  value =  vNew;
}

double getOptionValue(const vector<double> &sHat,
		      const vector<vector<double>>& vNew,
		      double S0, // initial stock value
		      double t0, // initial time
		      double T, // contract maturity date
		      double sigma, // volatility of default rate
		      double r, // risk free interest rate
		      double d, // risk free interest rate
		      double X, // risk free interest rate
		      int n
)
{
  // finish looping through time levels
  // return the interpolated value at lambda0
  double value=0.;
  double Shat0 = (S0*exp(-d*(T-t0))-X*exp(-r*(T-t0)))/sigma;
  // 
  if(Shat0 > sHat[sHat.size()-1] || Shat0 < sHat[0])
  {
    return outerSolution(S0,t0,T,X,r,d);
  }
  for(int ni=0;ni<=n;ni++)
  {
    value+=pow(sigma,ni+1)*interpolate(Shat0,sHat[1]-sHat[0],sHat,vNew[ni]);
  }
  return value;
} 

double getOptionValue(const vector<double> &sHat,
		      const vector<vector<double>>& vNew,
		      double S0, // initial stock value
		      double t0, // initial time
		      double T, // contract maturity date
		      double sigma, // volatility of default rate
		      double r, // risk free interest rate
		      double d, // risk free interest rate
		      double X // risk free interest rate
)
{
  return getOptionValue(sHat,vNew,S0,t0,T,sigma,r,d,X,vNew.size()-1);
} 

/* Template code for the Crank Nicolson Finite Difference
 */
int main()
{
  // parameters for the problem
  double t0=0.,T=0.5,sigma=0.2,r=0.05,d=0.15,X=2.;
  // n is the maximum number of terms to calculate
  // iMax is the number of time steps
  // jMax is the number of space steps
  // sMax is related to the domain of the grid
  int n=2,iMax=500,jMax=500,sMax=5;
  
  vector<double> sHat;
  vector<vector<double>> value;
  
  // calculate all the solutions C_i for 0<=i<=n
  caluculateOptionPrice(sHat,value,t0,T,r,X,n,iMax,jMax,sMax);
  
  vector<double> errors(n+1,0.);
  
  // print out solutions at different values of S
  ofstream output("test.csv");
  output.precision(12);
  for(int i=-100;i<=100;i++)
  {
    double S0 = X*exp(-r*(T-t0))*(1+i/100.);
    double analyticValue = putOptionPrice(S0,t0,X,r,d,sigma,T);
    output << S0;
    for(int ni=0;ni<=n;ni++)
    {
      double C_ni_Solution = getOptionValue(sHat,value,S0,t0,T,sigma,r,d,X,ni);
      output << " , " << C_ni_Solution;
      errors[ni] = max(errors[ni],fabs(C_ni_Solution-analyticValue));
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
