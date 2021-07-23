#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
/** \brief bisection algorithm to find the min of a function
 * \details See findMaxBisection for details.
 * @param f the value function, f(x) which is to be minimised.
 * @param xMin minimum value of x, which is lower bound of interval, \f$ a \f$.
 * @param xMax maximum value of x, which is upper bound of interval, \f$ b \f$.
 * @param tol user supplied tolerance, such that \f$ x^m \in [ x^* - tol(b-a) , x^* + tol(b-a) ] \f$ where  \f$ x^m \f$ is the estimate and  \f$ x^* \f$ is the true minimum.
 * @return returns the root and min value
 */
template<class T,class F>
static std::pair<T,T> goldenSearch(const F &f,T xMin,T xMax,T tol)
{
    static const T phi = (3-sqrt(5))/2.;
    // max number of iterations is related to the tolerance through the relation
    const int maxIter = std::max(int(log( tol )/log( 1.-phi )),1);
    // first get value at the left boundary
    // and put in as current min
    T xRoot=xMin; 
    T minValue=f(xRoot);
    //           // next check the +ve golden ratio value
    {
        T x = xMin + phi*(xMax-xMin);
        T temp=f(x);
        if(temp < minValue)
        {
            xRoot=x;
            minValue=temp;
        }
    }
    // the -ve one
    {
        T x = xMax - phi*(xMax-xMin);
        T temp=f(x);
        if(temp < minValue)
        {
            xRoot=x;
            minValue=temp;
        }
    }
    // and finally the right hand side boundary
    {
        T x = xMax;
        T temp=f(x);
        if(temp < minValue)
        {
            xRoot=x;
            minValue=temp;
        }
    }
    // now iterate using golden ratio
    for(int iter=1;iter<=maxIter;iter++)
    {
        // divide the interval [xMin,xMax] into 3 sections using dp
        T h=(xMax-xMin);
        // now evaluate the function at pMin+dp, pMin+2*dp, pMin+3*dp, given that we have already evaluated the function at xMin and xMax previously
        // -- check to see whether anyone is bigger than before
        T x;
        if(xRoot>xMin+0.5*h)
        {
            xMin=xMin+phi*h;
            x = xMax-phi*(xMax-xMin);
        }
        else
        {
            xMax=xMax-phi*h;
            x = xMin+phi*(xMax-xMin);
        }
        T temp=f(x);
        // if this is very close to min value then break
        if(temp < minValue)
        {
            xRoot=x;
            minValue=temp;
        }
        else if(std::fabs(temp-minValue)<tol)
            break;
        
    }
    return {xRoot,minValue};
}

// A generic lagrange interpolation function
double lagrangeInterpolation(const std::vector<double>& y,const std::vector<double>& x,double x0,unsigned int n)
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

struct GeneratorState
{
    // grid for the state variable
    std::vector<double> x;
    
    // value function on grid points
    std::vector<double> v;
    // value function on grid points
    std::vector<std::vector<double>> control;
    
    // value function operator
    double operator()(double x0){return lagrangeInterpolation(v,x,x0,2);}
    
    // dynamics of the state
    double f(double xt,double c) const {
        return std::max(cMin(xt),std::min(c,cMax(xt)));
    } 
    std::function<double(double)> cMin;
    std::function<double(double)> cMax;
    std::function<int(double)> eta;
    
    // cost function
    std::function<double(double,double,double)> Gamma;
    
};

struct Generator
{
    
    std::function<double(double)> EP;
    
    double tol;
    
    std::vector<GeneratorState> state;
    std::vector<double>  t;
    
    int solve();
    
    int outputPath(double x0,int u0,std::ostream& output=std::cout);
    
    template <class T>
    int setupState(int state,int n,int m,double tol,const T &params);
    
    GeneratorState& operator[](unsigned int U){return state[U];}
    
};


int main()
{
    int n=100;//number of steps per unit x,
    double T = 24;// number of hours
    int m=T*n;//number of timesteps,
    
    Generator G;
    G.tol=1.e-8;
    G.state.resize(2);
    
    // SETUP OFF STATE WITH INDEX 0
    // state 0 has points in interval [-1,0]
    G[0].x.resize(n+1,0.);
    for(int j=0;j<=n;j++)G[0].x[j] = -1. + (double)j/(double)n;// GRID
    G[0].v.resize(n+1,0.);
    // initialise value function
    G[0].control.clear();
    G[0].control.resize(n+1,std::vector<double>(m+1,0.));
    // define dynamics or minimum df/dx, normally negative number here
    G[0].cMin = [](double){return -1;};
    // define dynamics or maximum df/dx, normally positive number here
    G[0].cMax = [](double){return 1;};
    G[0].eta = [](double xtplus){
        // if xt would be positive, generator turned on
        if(xtplus>0.)return 1;  // 
        // otherwise stays off
        else return 0;
    };
    // cost function
    G[0].Gamma = [](double,double,double){return 0.;};
    
    // SETUP ON STATE WITH INDEX 1
    // state 0 has points in interval [0,1]
    G[1].x.resize(n+1,0.);
    for(int j=0;j<=n;j++)G[1].x[j] = (double)j/(double)n;// GRID
    G[1].v.resize(n+1,0.);
    // initialise value function
    G[1].control.clear();
    G[1].control.resize(n+1,std::vector<double>(m+1,0.));
    // define dynamics, any negative number here
    G[1].cMin = [](double){return -1;};
    // define dynamics, any positive number here
    G[1].cMax = [](double){return 1;};
    G[1].eta = [](double xtplus){
        // if xt would be negative, generator turns off
        if(xtplus<0.)return 0;
        // otherwise stays on
        else return 1;
    };
    // cost function
    G[1].Gamma = [](double xt,double pt,double c){
        return xt*pt - 30.;
    };
    
    
    G.t.resize(m+1);
    for(int k=0;k<=m;k++)
    {
        G.t[k] = k*T/m;
    }
    
    // setup electricity price function
    G.EP = [](double t){
        if(t<8.)return 20.;
        if(t<16.)return 40.;
        if(t<20.)return 60.;
        return 20.;
    };
    
    // run solver
    G.solve();
    
    //output path
    G.outputPath(0.,0);
    return 0;
}




int Generator::solve()
{
    
    // assign terminal conditions
    for(unsigned int i=0;i<state.size();i++)
    {
        for(unsigned int j=0;j<state[i].v.size();j++)
        {
            state[i].v[j]=0.;
        }
        
    }
    
    // timestep through solution
    for(int k=t.size()-2;k>=0;k--)
    {
        //             std::cout << "#####\n## Solve at time "<<t[k]<<"\n#####\n"<< std::endl;
        double dt = t[k+1]-t[k];
        
        std::vector<GeneratorState> Gold=state;
        
        for(unsigned int U=0;U<state.size();U++)
        {
            for(unsigned int j=0;j<state[U].x.size();j++)
            {
                double cMinStar=state[U].cMin(state[U].x[j]),cMaxStar=state[U].cMax(state[U].x[j]);
                auto objective = [&](double c){// calculate position of the characteristic using
                    double xHalfStar = state[U].x[j] + 0.5*state[U].f(state[U].x[j],c)*dt;
                    double xStar = state[U].x[j] + state[U].f(xHalfStar,c)*dt;
                    int uStar = state[U].eta(xStar);
                    
                    double temp=0.;
                    if( xStar < ( state[uStar].x.front() - tol ) )
                        temp-=fabs(xStar-state[uStar].x.front())/tol;
                    if( xStar > ( state[uStar].x.back() + tol ) )
                        temp-=fabs(xStar-state[uStar].x.back())/tol;
                    
                    // solving DV/Dt = g(x,t,c)
                    // so the simple explicit scheme is given by
                    // ( v(xStar,t+dt) - v(x,t) )/dt + g = 0
                    temp += Gold[uStar](xStar);
                    temp += state[U].Gamma(xHalfStar,EP(t[k]+0.5*dt),c)*dt;
                    return -temp;
                };
                std::pair<double,double> root = goldenSearch(objective,cMinStar,cMaxStar,tol);
                state[U].v[j] = -root.second;
                state[U].control[j][k]=root.first;
            }
        }
        
        if(k==0)
        {
            for(unsigned int U=0;U<state.size();U++)
            {
                int m=state[U].v.size()-1;
                for(unsigned int j=0;j<state[U].v.size();j+=std::max(m,1))
                {
                    std::cout << k << " " << U << " " << state[U].x[j] << " " <<  state[U].v[j]  << " " << state[U].control[j][k] <<std::endl;
                }
            }
        }
    }
    return 0;
};

int Generator::outputPath(double x0,int u0,std::ostream& output)
{
    double xt=x0;
    double ut=u0;
    double vt=0.;
    
    for(unsigned int k=0;k<t.size()-1;k++)
    {
        double dt=t[k+1]-t[k];
        // get control from nearest integer
        double dx;
        if(state[ut].x.size()>1)
            dx=state[ut].x[1]-state[ut].x[0];
        else
            dx=1.;
        int jStar = int((xt - state[ut].x[0])/dx);
        
        jStar = std::min((unsigned long)(jStar),state[ut].x.size()-1);
        int jStarPlus = std::min((unsigned long)(jStar+1),state[ut].x.size()-1);
        double ct = state[ut].control[jStar][k];
        
        if(jStar!=jStarPlus){
            ct+=(xt-state[ut].x[jStar])*(state[ut].control[jStarPlus][k]-ct)/dx;
        }
        
        double xHalfStar = xt + 0.5*state[ut].f(xt,ct)*dt;
        vt = vt + state[ut].Gamma(xHalfStar,EP(t[k]+0.5*dt),ct)*dt;
        xt = xt + state[ut].f(xHalfStar,ct)*dt;        
        ut = state[ut].eta(xt);
        
        xt = std::max(state[ut].x.front(),std::min(xt,state[ut].x.back()));
        
        output << t[k+1] <<" "<< ut <<" "<< xt << " "<< vt << " "<< ct << std::endl;
    }
    return 0;
}
