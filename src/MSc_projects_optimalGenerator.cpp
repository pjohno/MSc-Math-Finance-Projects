#include "MSc_projects_optimalGenerator.hpp"
#include "math60082_markupTable.hpp"
#include "math60082_mMinima.hpp"
using namespace std;

namespace MSC_PROJECTS
{
    
    double GeneratorState::f(double xt,double c) const {
        return std::max(cMin(xt),std::min(c,cMax(xt)));
    } 
    
    int Generator::defaultSetup(int n,int m,double T,double tol_){
        
        state.clear();
        // create a 4 state model
        state.resize(4);
        tol=tol_;
        
        // state 0 represents "off" state 
        // x goes from -1 to 0
        // -1 < c < 1
        // no cost
        state[0].x.resize(n+1,0.);
        for(int i=0;i<=n;i++)state[0].x[i] = -1. + (double)i/(double)n;// initialise value function
        state[0].v.resize(n+1,0.);
        // initialise value function
        state[0].control.resize(n+1,std::vector<double>(m+1,0.));
        // define dynamics
        state[0].cMin = [](double){return -1;};
        // define dynamics
        state[0].cMax = [](double){return 1.;};
        // start warming up
        state[0].eta = [&](double xtplus,double ct){
            // if xt would be positive start warming up
            if(xtplus>-tol)return 1;
            else return 0;
        };
        // zero cost function
        state[0].Gamma = [](double xt,double pt,double c){return 0.;};
        
        // state 1 represents warming up
        // 0 <= x <= 1 combined with c=1 means that
        // it takes 1 unit of time to move through
        state[1].x.resize(n+1);
        for(int i=0;i<=n;i++)state[1].x[i] = (double)i/(double)n;
        // initialise value function
        state[1].v.resize(n+1,0.);
        // initialise value function
        state[1].control.resize(n+1,std::vector<double>(m+1,0.));
        // define dynamics
        state[1].cMin = [](double){return 1.;};
        // define dynamics
        state[1].cMax = [](double){return 1.;};
        state[1].eta = [&](double xtplus,double ct){
            // switch to "on" state 2 if reach x=1
            if(xtplus>1.-tol)return 2;
            else return 1;
        };
        // cost function
        state[1].Gamma = [](double xt,double pt,double c){
            // if warming up, some fixed costs incurred
            return -5;
        };
        
        // state 2 represents warming up
        // 1 <= x <= 2 combined with c_max=1 means that
        // it takes 1 unit of time to move through from 
        // maximum power to minimum. You could adjust cMin/cMax
        // to include ramp rates as they can be function of "x"
        state[2].x.resize(n+1);
        for(int i=0;i<=n;i++)state[2].x[i] = 1.+(double)i/(double)n;
        // initialise value function
        state[2].v.resize(n+1,0.);
        // initialise value function
        state[2].control.resize(n+1,std::vector<double>(m+1,0.));
        // define dynamics
        state[2].cMin = [](double xt){return -1.;};
        // define dynamics
        state[2].cMax = [](double xt){return 1.;};
        state[2].eta = [&](double xtplus,double ct){
            // switch to "warming down" state 3 if reach x=1
            if(xtplus<1.+tol)return 3;
            else return 2;
        };
        // cost function
        state[2].Gamma = [](double xt,double pt,double c){
            // 
            return (xt-1.)*pt - 30.;
        };
        
        // state 1 represents warming up
        // 0 <= x <= 1 combined with c=1 means that
        // it takes 1 unit of time to move through
        state[3].x.resize(n+1);
        for(int i=0;i<=n;i++)state[3].x[i] = (double)i/(double)n;
        // initialise value function
        state[3].v.resize(n+1,0.);
        // initialise value function
        state[3].control.resize(n+1,std::vector<double>(m+1,0.));
        // define dynamics
        state[3].cMin = [](double){return -1.;};
        // define dynamics
        state[3].cMax = [](double){return -1.;};
        state[3].eta = [&](double xtplus,double ct){
            // switch to "off" state 0 if reach x=0
            if(xtplus<tol)return 0;
            else return 3;
        };
        // cost function
        state[3].Gamma = [](double xt,double pt,double c){
            // if warming down, some fixed costs incurred
            return -5;
        };
        
        t.resize(m+1);
        for(int k=0;k<=m;k++)
        {
            t[k] = k*T/m;
        }
        
        return 0;
    }
    
    int Generator::solve()
    {
        using MATH60082::goldenSearch;
        
        // assign terminal conditions
        for(unsigned int U=0;U<state.size();U++)
        {
            for(unsigned int j=0;j<state[U].v.size();j++)
            {
                state[U].v[j]=0.;
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
                        int uStar = state[U].eta(xStar,c);
                        
                        double temp=0.;
                        if( xStar < ( state[uStar].x.front() - tol ) )
                            temp-=fabs(xStar-state[uStar].x.front())/tol;
                        if( xStar > ( state[uStar].x.back() + tol ) )
                            temp-=fabs(xStar-state[uStar].x.back())/tol;
                        
                        // solving DV/Dt = g(x,t,c)
                        // so the simple explicit scheme is given by
                        // ( v(xStar,t+dt) - v(x,t) )/dt + g = 0
                        temp += Gold[uStar](xStar);
                        temp += state[U].Gamma(xHalfStar,electricityPrice(t[k]+0.5*dt),c)*dt;
                        return -temp;
                    };
                    std::pair<double,double> root = goldenSearch(objective,cMinStar,cMaxStar,tol);
                    state[U].v[j] = -root.second;
                    state[U].control[j][k]=root.first;
                }
            }
            
            if(k==0)
            {
                MATH60082::tableRow( "k" , "U_t" , "x_t" ,"v_t" , "c_t");
                MATH60082::emptyTableRow( 5 );
                for(unsigned int U=0;U<state.size();U++)
                {
                    int m=state[U].v.size()-1;
                    for(unsigned int j=0;j<state[U].v.size();j+=std::max(m,1))
                    {
                        MATH60082::tableRow( k ,  U , state[U].x[j] , state[U].v[j] , state[U].control[j][k] );
                    }
                }
            }
        }
        return 0;
    };
    
    int Generator::outputPath(double x0,int u0,std::ostream& output,bool toMarkup)
    {
        double xt=x0;
        double ut=u0;
        double vt=0.;
        if(toMarkup){
        MATH60082::tableRow( "t_k" , "U_t" , "x_t" ,"v_t" , "c_t", "p_t");
        MATH60082::emptyTableRow( 6 );
        MATH60082::tableRow( t[0] , ut , xt ,vt , 0. , electricityPrice(t[0]));
        }
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
            vt = vt + state[ut].Gamma(xHalfStar,electricityPrice(t[k]+0.5*dt),ct)*dt;
            xt = xt + state[ut].f(xHalfStar,ct)*dt;        
            ut = state[ut].eta(xt,ct);
            
            xt = std::max(state[ut].x.front(),std::min(xt,state[ut].x.back()));
            if(toMarkup)
                MATH60082::tableRow( t[k+1] , ut , xt ,vt , ct, electricityPrice(t[0]));
            else
                output << t[k+1] <<" "<< ut <<" "<< xt << " "<< vt << " "<< ct << " " << electricityPrice(t[k+1])<< endl;
        }
        return 0;
    }
    
}

