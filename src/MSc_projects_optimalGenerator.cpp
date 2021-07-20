#include "MSc_projects_optimalGenerator.hpp"
#include "math60082_markupTable.hpp"
#include <boost/math/tools/minima.hpp>
using namespace std;

namespace MSC_PROJECTS
{
    
    double GeneratorState::f(double xt,double c) const {
        return std::max(cMin(xt),std::min(c,cMax(xt)));
    } 
    
    template <>
    int Generator::setupState<GeneratorDefaultOff>(int state,int n,int m,double tol,const GeneratorDefaultOff &p)
    {
        // state 0 has only one point
        G[state].x = Vector(1,0.);
        // initialise value function
        G[state].v = Vector(1,0.);
        // initialise value function
        G[state].control = Matrix(1,m+1,0.);
        // define dynamics
        G[state].cMin = [p](double){
            // nowhere to go, just return 0
            return p.cMin;
        };
        // define dynamics
        G[state].cMax = [p](double){
            // nowhere to go, just return 0
            return p.cMax;
        };
        G[state].eta = [tol,p](double xtplus){
            // if xt would be positive start warming up
            if(xtplus>tol)return 1;
            else return 0;
        };
        // cost function
        G[state].Gamma = [p](double xt,double pt,double c){
            // if turned off, no need to include costs
            return p.cost;
        };
        return 0;
    }
  
    template <>
    int Generator::setupState<GeneratorDefaultSyncing>(int state,int n,int m,double tol,const GeneratorDefaultSyncing &p)
    {
        // state 0 has only one point
        G[state].x = Vector(n+1);
        for(int i=0;i<=n;i++)
        {G[state].x[i] = p.xMin + i*(p.xMax-p.xMin)/n;}
        // initialise value function
        G[state].v = Vector(n+1,0.);
        // initialise value function
        G[state].control = Matrix(n+1,m+1,0.);
        // define dynamics
        G[state].cMin = [p](double){
            // rate of 1
            return p.cMin;
        };
        // define dynamics
        G[state].cMax = [p](double){
            // rate of 1
            return p.cMax;
        };
        G[state].eta = [tol,p](double xtplus){
            // while xt 
            if(xtplus<p.xMax-tol)return 1;
            else return 2;
        };
        // cost function
        G[state].Gamma = [p](double xt,double pt,double c){
            // if warming up, some fixed costs incurred
            return p.cost;
        };
        return 0;
    }
   
    template <>
    int Generator::setupState<GeneratorDefaultOn>(int state,int n,int m,double tol,const GeneratorDefaultOn &p)
    {
        // state 3 has only one point
        G[state].x = Vector(n+1);
        for(int i=0;i<=n;i++)
        {G[state].x[i] = p.xMin + i*(p.xMax-p.xMin)/n;}
// initialise value function
        G[state].v = Vector(n+1,0.);
        // initialise value function
        G[state].control = Matrix(n+1,m+1,0.);
        // define dynamics
        G[state].cMin = [p](double){
            // rate of 1
            return p.cMin;
        };
        // define dynamics
        G[state].cMax = [p](double){
            // rate of 1
            return p.cMax;
        };
        G[state].eta = [tol,p](double xtplus){
            // while xt 
            if(xtplus>p.xMin+tol)return 2;
            else return 3;
        };
        // cost function
        G[state].Gamma = p.cost;
//         [](double xt,double pt,double c){
//             // if warming up, some fixed costs incurred
//             return pt;
//         };
        
        return 0;
    }
    
    template <>
    int Generator::setupState<GeneratorDefaultDeSyncing>(int state,int n,int m,double tol,const GeneratorDefaultDeSyncing &p)
    {
        // state 3 has only one point
        G[state].x = Vector(n+1);
        for(int i=0;i<=n;i++)
        {G[state].x[i] = p.xMin + i*(p.xMax-p.xMin)/n;}
        // initialise value function
        G[state].v = Vector(n+1,0.);
        // initialise value function
        G[state].control = Matrix(n+1,m+1,0.);
        // define dynamics
        G[state].cMin = [p](double){
            // rate of 1
            return p.cMin;
        };
        G[state].cMax = [p](double){
            // rate of 1
            return p.cMax;
        };
        G[state].eta = [tol,p](double xtplus){
            // while xt 
            if(xtplus>p.xMin+tol)return 3;
            else return 0;
        };
        // cost function
        G[state].Gamma = [p](double xt,double pt,double c){
            // if warming up, some fixed costs incurred
            return p.cost;
        };
        
        return 0;
    }
    
    
    Generator::Generator(int n,int m,double T,double tol_){
        
        G.clear();
        G.resize(4);
        tol=tol_;
        
        setupState<GeneratorDefaultOff>(0,n,m,tol,{0.,0.,0.,1.,0.});
        setupState<GeneratorDefaultSyncing>(1,n,m,tol,{0.,1.,0.,1.,-1.});
        setupState<GeneratorDefaultOn>(2,n,m,tol,{1.,2.,-1.,1.,[](double xt,double pt,double c){
            // if warming up, some fixed costs incurred
            return xt*pt;
        }});
        setupState<GeneratorDefaultDeSyncing>(3,n,m,tol,{0.,1.,-1.,0.,-1.});
        
        t = Vector(m+1);
        for(int i=0;i<=m;i++)
        {t[i] = i*T/m;}
        
    }
    
    int Generator::solve()
    {
        using boost::math::tools::brent_find_minima;
        
        // assign terminal conditions
        for(unsigned int i=0;i<G.size();i++)
        {
            for(unsigned int j=0;j<G[i].v.size();j++)
            {
                G[i].v[j]=0.;
            }
            
        }
        
        // timestep through solution
        for(int k=t.size()-2;k>=0;k--)
        {
            //             std::cout << "#####\n## Solve at time "<<t[k]<<"\n#####\n"<< std::endl;
            double dt = t[k+1]-t[k];
            
            std::vector<GeneratorState> Gold=G;
            
            for(unsigned int i=0;i<G.size();i++)
            {
                for(unsigned int j=0;j<G[i].v.size();j++)
                {
                    double uMinStar=G[i].cMin(G[i].x[j]),uMaxStar=G[i].cMax(G[i].x[j]);
                    std::pair<double,double> root = brent_find_minima([&](double u){
                                                                            // calculate position of the characteristic using
                                                                            double xHalfStar = G[i].x[j] + 0.5*G[i].f(G[i].x[i],u)*dt;
                                                                            double xStar = G[i].x[j] + G[i].f(xHalfStar,u)*dt;
                                                                            int iStar = G[i].eta(xStar);
                                                                            
                                                                            double temp=0.;
                                                                            if( xStar<G[iStar].x.front()-tol )
                                                                                temp-=fabs(xStar-G[iStar].x.front())/tol;
//                                                                                 temp-=1./tol;
                                                                            if( xStar>G[iStar].x.back()+tol )
                                                                                temp-=fabs(xStar-G[iStar].x.back())/tol;
//                                                                                 temp-=1./tol;
                                                                        
                                                                        // solving DV/Dt = g(x,t,u)
                                                                        // so the simple explicit scheme is given by
                                                                        // ( v(xStar,t+dt) - v(x,t) )/dt + g = 0
                                                                        temp += Gold[iStar](xStar);
                                                                        temp += G[i].Gamma(xHalfStar,EP(t[k]+0.5*dt),u)*dt;
                                                                        return -temp;
                                                                        },
                                                                        uMinStar,uMaxStar,20);
                    G[i].v[j] = -root.second;
                    G[i].control[j][k]=root.first;
                }
            }
            
            if(k==0)
            {
                for(unsigned int i=0;i<G.size();i++)
                {
                    int m=G[i].v.size()-1;
                    for(unsigned int j=0;j<G[i].v.size();j+=std::max(m,1))
                    {
                        double xHalfStar = G[i].x[j] + 0.5*G[i].f(G[i].x[i],G[i].control[j][k])*dt;
                        double xStar = G[i].x[j] + G[i].f(xHalfStar,G[i].control[j][k])*dt;
                        std::cout << k << " " <<  i << " " << G[i].x[j] << " " << G[i].v[j] << " " << G[i].control[j] <<  " " << G[i].eta(xStar) << std::endl;
                    }
                    std::cout << std::endl;
                }
            }
        }
        return 0;
    };
    
    int Generator::outputPath(double x0,int u0)
    {
        double xt=x0;
        double ut=u0;
        double vt=0.;
        MATH60082::tableRow( t[0] , 0. , xt ,ut , vt);
        for(unsigned int k=0;k<t.size()-1;k++)
        {
            double dt=t[k+1]-t[k];
            // get control from nearest integer
            double dx;
            if(G[ut].x.size()>1)
                dx=G[ut].x[1]-G[ut].x[0];
            else
                dx=1.;
            int jStar = int((xt - G[ut].x[0])/dx);
            jStar = std::min((unsigned int)(jStar),G[ut].x.size()-1);
            int jStarPlus = std::min((unsigned int)(jStar+1),G[ut].x.size()-1);
            double ct = G[ut].control[jStar][k];
            
            if(jStar!=jStarPlus){
                ct+=(xt-G[ut].x[jStar])*(G[ut].control[jStarPlus][k]-ct)/dx;
                cout << jStar << " " << jStarPlus << " " << G[ut].control[jStar][k] << " " << ct << " " << G[ut].control[jStarPlus][k] << endl;
            }
            
            
            
            double xHalfStar = xt + 0.5*G[ut].f(xt,ct)*dt;
            vt = vt + G[ut].Gamma(xHalfStar,EP(t[k]+0.5*dt),ct)*dt;
            xt = xt + G[ut].f(xHalfStar,ct)*dt;        
            ut = G[ut].eta(xt);
            
            xt = std::max(G[ut].x.front(),std::min(xt,G[ut].x.back()));
            MATH60082::tableRow( t[k+1] , ct , xt ,ut , vt);
        }
        return 0;
    }
    
}

