#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
using namespace std;

class MVector
{
    unsigned int N;
    double* v;
public:
    explicit MVector():N(0),v(nullptr){}
    explicit MVector(unsigned int n):N(n),v(new double[n]){}
    explicit MVector(unsigned int n,double x):N(n),v(new double[n]){for(unsigned int i=0;i<N;i++)v[i]=x;}
    ~MVector(){if(v!=nullptr) delete [] v;}
    // normal copies
    MVector(const MVector& X):N(X.N),v(new double[X.N]){std::copy(X.v,X.v+X.N,v);}
    MVector& operator=(const MVector &X){
        if(v!=nullptr && N!=X.N){delete [] v;}
        if(N!=X.N){v = new double [X.N];}
        N = X.N;
        std::copy(X.v,X.v+X.N,v);
        return *this;
    }
    // c++11 initialiser and move construct, move=
    explicit MVector(std::initializer_list<double> list):N(list.size()),v(new double[list.size()]){auto  vPos=v;for(auto lPos=list.begin();lPos!=list.end();lPos++,vPos++)*vPos=*lPos;}
    MVector(MVector&& X):N(X.N),v(X.v){X.N=0;X.v=nullptr;}
    MVector& operator=(MVector&& X){
        if(v!=nullptr){
            if(v==X.v){return *this;}
            delete [] v;
        }
        N=X.N;v=X.v;X.N=0;X.v=nullptr;
        return *this;
    }
    
    double& operator[](int index){return v[index];}
    double operator[](int index) const {return v[index];}
    unsigned int size() const {return N;}
    // resize vector
    void resize(int n){
        if(v!=nullptr) delete [] v;
        v=new double[n];
        N=n;
    }
    void resize(int n,double x){
        resize(n);
        double *vPos=v;
        for(unsigned int i=0;i<N;i++)*vPos++=x;
    }
    // return array
    double* returnArray(){return v;}
    const double* returnArray() const {return v;}
};

std::ostream& operator<<(std::ostream &output,const MVector &X)
{
  //output << "( ";
  unsigned int morethanone=std::min((unsigned int)(1),X.size());
  for(unsigned int i=0;i<morethanone;i++)
    output << X[i];
  for(unsigned int i=morethanone;i<X.size();i++)
    output << " , " << X[i];
  //output << " )";
  return output;
}

MVector operator*(const double& lhs,const MVector &rhs)
{
    MVector temp(rhs);
    for(unsigned int i=0;i<rhs.size();i++)
        temp[i]*=lhs;
    return temp;
}

MVector operator*(const double& lhs,MVector&& rhs)
{
    MVector temp(std::move(rhs));
    double *tempPos=temp.returnArray();
    for(unsigned int i=0;i<temp.size();i++)
        (*tempPos++)*=lhs;
    return temp;
}

MVector operator*(const MVector& lhs,const double &rhs)
{
    MVector temp(lhs);
    for(unsigned int i=0;i<lhs.size();i++)
        temp[i]*=rhs;
    return temp;
}

MVector operator/(const MVector& lhs,const double &rhs)
{
    MVector temp(lhs);
    for(unsigned int i=0;i<lhs.size();i++)
        temp[i]/=rhs;
    return temp;
}

MVector operator+(const MVector& lhs,const MVector &rhs)
{
    MVector temp(lhs);
    double *tempPos=temp.returnArray();
    const double *rPos=rhs.returnArray();
    for(unsigned int i=0;i<lhs.size();i++)
        (*tempPos++)+=(*rPos++);
    return temp;
}

MVector operator+(const MVector& lhs,MVector&& rhs)
{
    MVector temp(std::move(rhs));
    double *tempPos=temp.returnArray();
    const double *lPos=lhs.returnArray();
    for(unsigned int i=0;i<lhs.size();i++)
        (*tempPos++)+=(*lPos++);
    return temp;
}

MVector operator-(const MVector& lhs,const MVector &rhs)
{
  MVector temp(lhs);
  for(unsigned int i=0;i<lhs.size();i++)
    temp[i]-=rhs[i];
  return temp;
}

MVector operator-(const MVector& lhs,MVector&& rhs)
{
  MVector temp(std::move(rhs));
  double *tempPos=temp.returnArray();
  const double *lPos=lhs.returnArray();
  for(unsigned int i=0;i<lhs.size();i++)
  {
    *tempPos*=-1;
    (*tempPos++)+=(*lPos++);
  }
  return temp;
}

// n is number of steps
// t_0 = a
// t_n = b
// t_i = a + ih where h=(b-a)/n 
// x_0 = alpha
// this will solve to find x_n
template<class SCALAR,class VECTOR,class FUNCTION>
VECTOR RK4MethodTemplate(int n,SCALAR a,SCALAR b,
			    VECTOR alpha,const FUNCTION &f)
{
  // local variables
  SCALAR h,t;
  VECTOR x;
  // intialise values
  h=(b-a)/(SCALAR)(n);
  t=a;
  x=alpha;
  // implement Euler's method
  for(int i=0;i<n;i++)
  {
    t = a + i*h; // update value of t to t_i
    VECTOR k1 = h*f(x,t); // update x to x_{i+1}
    VECTOR k2 = h*f(x+0.5*k1,t+0.5*h); // update x to x_{i+1}
    VECTOR k3 = h*f(x+0.5*k2,t+0.5*h); // update x to x_{i+1}
    VECTOR k4 = h*f(x+k3,t+h); // update x to x_{i+1}
    x = x + 1./6.*(k1+k4+(2.*(k2+k3))); 
  }
  return x; // returns x_n
}

int main()
{
    
    // parameters from Gueant et al 2012
    // figure 1
    int qMax=6;
    double mu=0.0;
    double sigma=0.3;
    double A = 0.1;
    double k=0.3;
    double gamma = 0.05;
    double b=3;
    double T= 300;
    
    // number of time observations
    int n=100;
    // time step
    double dT = T/n;
    
    // store the value of the omega and delta at each time step and q value
    vector<MVector> omega(n+1,MVector(qMax+1)),delta(n+1,MVector(qMax+1));
    
    // initialise the solution at t=T
    for(int q=0;q<=qMax;q++)
    {
        // from initial condition
        omega[n][q] = exp(-k*q*b);
        // from formula for delta
        if(q>0)
            delta[n][q] = 1/k*log(omega[n][q]/omega[n][q-1]) + 1/gamma*log(1+gamma/k);
        else
            delta[n][q] = 0.;
    }
    
    // now use a numerical integration to find the value at T_i given the value at T_{i+1}
    for(int i = n-1 ; i>=0 ; i--)
    {
        // some constants in the equation
        double alpha = k/2.*gamma*sigma*sigma;
        double beta = k*mu;
        double eta = A*pow(1+gamma/k,-(1+k/gamma));
        // now solve 
        //   dw(q,t)/dt = (alpha q^2 - beta q) w(q,t) - eta w(q-1,t)
        // with initial condition 
        //   w(q,T_i) = omega(q,T_i)
        // so that
        //   omega(q,T_{i-1}) = w(q,T_{i-1})
        omega[i] = RK4MethodTemplate(100,i*dT,(i-1)*dT,omega[i+1],
                                     [&]
                                     (const  MVector &w,double t)
                                     {
                                         MVector F(qMax+1);
                                         F[0] = 0.;
                                         for(int q=1;q<=qMax;q++)
                                             F[q] = (alpha*q*q - beta*q)*w[q] - eta*w[q-1];
                                         return F;
                                     }
        );
        
        // We can then calculate the value of the optimal ask price
        for(int q=1;q<=qMax;q++)
            delta[i][q] = 1/k*log(omega[i][q]/omega[i][q-1]) + 1/gamma*log(1+gamma/k);
        delta[i][0] = 0.;
    }
    
    // and output results to file
    ofstream output("test.csv");
    for(int i = 0 ; i<=n ; i++)
    {
        output << i*dT << " , " << omega[i] << " , " << delta[i] << endl;
    }
    
}

