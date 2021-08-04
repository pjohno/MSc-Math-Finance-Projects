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
    
    int Generator::karolisExampleSetup(int stepsPerHour,double T,double tol_){
        // get dt
        double dt=1./stepsPerHour;
        // total timestep m= stepsperhour*T
        int m=stepsPerHour*T;
        t.resize(m+1);
        for(int k=0;k<=m;k++)
        {
            t[k] = k*dt;
        }
        
        state.clear();
        // create a 6 state model
        state.resize(6);
        tol=tol_;
        
        // NDZ is 12 hours or 720 mins
        double NDZ=12.;
        // NDZ is 6 hours or 360 mins
        double MZT=6.;
        // l is labour costs -- set to zero
        double l=0;
        // efficiency
        double eta_min=0.2;
        // efficiency
        double eta_max=0.5;
        // efficiency
        double gamma=0.413185;
  
        // export limits
        double SEL=50.;
        // export limits
        double MEL=90.;
        
        // fuel cost
        double lambda=25.;
        // zi is fuel cost running at full power
        double zi=lambda*MEL/eta_max;
        // percentage additional syncing cost
        double psi_fcOM=0.2;
        // percentage fixed cost
        double psi_fc=0.2;
 
        // run up rates
        double RURE1=1.*60.;
        // run up rates
        double RURE2=0.1*60.;
        // run up rates
        double RURE3=2.*60.;
        // run up limits
        double RUEE2=30.;
        // run up limits
        double RUEE3=60.;
       
        // run down rates
        double RDRE1=2.*60.;
        // run up rates
        double RDRE2=0.1*60.;
        // run up rates
        double RDRE3=1.*60.;
        // run up limits
        double RDEE2=60.;
        // run up limits
        double RDEE3=30.;
        
        {
            // state x<-NDZ represents "off" state 
            // x goes from -NDZ-1 to -NDZ
            // -1 < c < 1
            // no cost
            int n=5;
            
            int current=0;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            
            state[current].x.resize(n+1,0.);
            for(int j=0;j<=n;j++)state[current].x[j] = -NDZ-1. + (double)j/(double)n;// initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [](double){return -1;};
            // define dynamics
            state[current].cMax = [](double){return 1.;};
            // start warming up
            state[current].eta = [this,NDZ](double xtplus,double ct){
                // if xt would be positive start syncing up
                if(xtplus>-NDZ-tol)return 1;
                else return 0;
            };
            // just labour cost function
            state[current].Gamma = [l](double xt,double pt,double c){return -l;};
        }
        {
            int n=NDZ*stepsPerHour;
            double dx=dt;
            // state 1 represents syncing up
            // 0 <= x <= 1 combined with c=1 means that
            // it takes 12 units of time to move through
            int current=1;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            state[current].x.resize(n+1);
            for(int i=0;i<=n;i++)state[current].x[i] = -NDZ + i*dx;
            // initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [](double){return 1.;};
            // define dynamics
            state[current].cMax = [](double){return 1.;};
            state[current].eta = [this](double xtplus,double ct){
                // switch to "on" state 2 if reach x=0
                if(xtplus>-tol)return 2;
                else return 1;
            };
            // cost function
            state[current].Gamma = [zi,psi_fcOM,psi_fc,l](double xt,double pt,double c){
                // if warming up, some fixed costs incurred
                return -zi*(psi_fcOM+psi_fc)-l;
            };
        }
        
        {
            // state 2 represents warming up to the SEL
            // 0 <= x <= SEL combined with RURE1=1./60 means that
            // it takes 1.6666 unit of time to move through
            
            // get approx dx
            double dx=RURE1*dt;
            // get number of points in grid
            int n = SEL/dx + 0.5;
            dx = SEL/double(n);
            
            int current=2;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            state[current].x.resize(n+1);
            for(int i=0;i<=n;i++)state[current].x[i] = i*dx;
            // initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [RURE1,RURE2,RUEE2,SEL,dt](double xt){
                if(xt<RUEE2)
                    return RURE1;
                else
                {
                    if(xt+RURE2*dt>SEL)
                        return (xt-SEL)/dt;
                    else
                        return RURE2;
                }
            };
            // define dynamics
            state[current].cMax = [RURE1,RURE2,RUEE2,SEL,dt](double xt){
                if(xt<RUEE2)
                    return RURE1;
                else
                {
                    if(xt+RURE2*dt>SEL)
                        return (xt-SEL)/dt;
                    else
                        return RURE2;
                }
            };
            state[current].eta = [this,SEL](double xtplus,double ct){
                // switch to "on" state 3 if reach x=1
                if(xtplus>SEL-tol)return 3;
                else return 2;
            };
            // cost function
            state[current].Gamma = [this,zi,psi_fcOM,l,psi_fc,lambda,gamma,eta_min,eta_max,MEL](double xt,double pt,double c){
                if(xt<tol)return -zi*(psi_fcOM+psi_fc)-l;
                double percentageLoad=xt/MEL;
                double  efficiency=eta_max*pow(percentageLoad,gamma);
                double fuelInput;
                
                if(efficiency<eta_min)
                    fuelInput = xt/eta_min;
                else
                    fuelInput = MEL/eta_max*pow(percentageLoad,1-gamma);
                
                // if warming up, some fixed costs incurred
                return xt*pt-lambda*fuelInput-zi*psi_fc-l;
            };
        }
        
        {
            // state 3 represents on state between SEL and MEL
            // SEL <= x <= MEL combined with RURE2/3
            // get approx dx
            double dx=RURE3*dt;
            // get number of points in grid
            int n = (MEL-SEL)/dx + 0.5;
            dx = (MEL-SEL)/double(n);
            
            int current=3;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            state[current].x.resize(n+1);
            for(int i=0;i<=n;i++)state[current].x[i] = SEL+i*dx;
            // initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [RDRE1,RDRE2,RDRE3,RDEE2,RDEE3](double xt){
                if(xt>RDEE2)
                    return -RDRE1;
                else if(xt>RDEE3)
                    return -RDRE2;
                else
                    return -RDRE3;
            };
            // define dynamics
            state[current].cMax =[RURE1,RURE2,RURE3,RUEE2,RUEE3](double xt){
                if(xt<RUEE2)
                    return RURE1;
                else if(xt<RUEE3)
                    return RURE2;
                else
                    return RURE3;
            };
            state[current].eta = [this,SEL](double xtplus,double ct){
                // switch to "warming down" state 3 if reach x=SEL
                if(xtplus<SEL-tol)return 4;
                else return 3;
            };
            // cost function
            state[current].Gamma =  [this,psi_fcOM,zi,l,psi_fc,lambda,gamma,eta_min,eta_max,MEL](double xt,double pt,double c){
                if(xt<tol)return -zi*(psi_fcOM+psi_fc)-l;
                double percentageLoad=xt/MEL;
                double  efficiency=eta_max*pow(percentageLoad,gamma);
                double fuelInput;
                
                if(efficiency<eta_min)
                    fuelInput = xt/eta_min;
                else
                    fuelInput = MEL/eta_max*pow(percentageLoad,1-gamma);
                
                // if warming up, some fixed costs incurred
                return xt*pt-lambda*fuelInput-zi*psi_fc-l;
            };
        }
        
        {
            // state 4 represents on state between SEL and 0
            // 0 <= x <= SEL combined with RDRE2/3
            // get approx dx
            double dx=RDRE3*dt;
            // get number of points in grid
            int n = SEL/dx+0.5;
            dx = SEL/double(n);
            
            int current=4;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            state[current].x.resize(n+1);
            for(int i=0;i<=n;i++)state[current].x[i] = i*dx;
            // initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [RDRE1,RDRE2,RDRE3,RDEE2,RDEE3,dt](double xt){
                if(xt>RDEE2)
                    return -RDRE1;
                else if(xt>RDEE3)
                    return -RDRE2;
                else
                {
                    if(xt-RDRE3*dt<0.)
                        return -xt/dt;
                    else
                        return -RDRE3;
                }
            };
            // define dynamics
            state[current].cMax = [RDRE1,RDRE2,RDRE3,RDEE2,RDEE3,dt](double xt){
                if(xt>RDEE2)
                    return -RDRE1;
                else if(xt>RDEE3)
                    return -RDRE2;
                else
                {
                    if(xt-RDRE3*dt<0.)
                        return -xt/dt;
                    else
                        return -RDRE3;
                }
            };
            state[current].eta = [this](double xtplus,double ct){
                // switch to "syncing off" state 5 if reach x=0
                if(xtplus<tol)return 5;
                else return 4;
            };
            // cost function
            state[current].Gamma =  [this,zi,l,psi_fcOM,psi_fc,lambda,gamma,eta_min,eta_max,MEL](double xt,double pt,double c){
                if(xt<tol)return -zi*(psi_fcOM+psi_fc)-l;
                double percentageLoad=xt/MEL;
                double  efficiency=eta_max*pow(percentageLoad,gamma);
                double fuelInput;
                
                if(efficiency<eta_min)
                    fuelInput = xt/eta_min;
                else
                    fuelInput = MEL/eta_max*pow(percentageLoad,1-gamma);
                
                // if warming up, some fixed costs incurred
                return xt*pt-lambda*fuelInput-zi*psi_fc-l;
            };
        }
        
        {
            int n=NDZ*stepsPerHour;
            double dx=dt;
            // state 5 represents syncing off
            // -NDZ <= x <=  combined with c=-1 means that
            // it takes 12 units of time to move through
            int current=5;
            cout << " SETUP STATE " << current << " WITH n=" << n << endl;
            state[current].x.resize(n+1);
            for(int i=0;i<=n;i++)state[current].x[i] = -NDZ+i*dx;
            // initialise value function
            state[current].v.resize(n+1,0.);
            // initialise value function
            state[current].control.resize(n+1,std::vector<double>(m+1,0.));
            // define dynamics
            state[current].cMin = [](double xt){return -1.;};
            // define dynamics
            state[current].cMax = [MZT](double xt){
                if(xt<-MZT/2.)
                    return 1.;
                else 
                    return -1.;
            };
            state[current].eta = [this,NDZ,MZT](double xtplus,double ct){
                // switch to "on" state 2 if reach x=0
                if((xtplus<-MZT/2.) && ct>0.)return 1;
                if(xtplus<-NDZ+tol)return 0;
                else return 5;
            };
            // cost function
            state[current].Gamma = [zi,psi_fcOM,psi_fc,l](double xt,double pt,double c){
                // if warming up, some fixed costs incurred
                return -zi*(psi_fcOM+psi_fc)-l;
            };
        }
        
        electricityPrice = [](double t){
            if(t<24.)return 100.;
            if(t<28.)return -700.;
            if(t<48.)return 200.;
            if(t<80.)return 100.;
            if(t<85.)return 90.;
            if(t<90.)return 0.;
            return 80.;
        };
//          electricityPrice = [](double t){
//             if(t<60)
//                 return 100.;
//             else if(t<100)
//                 return -700.;
//             else 
//                 return 100.;
//         };
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
                MATH60082::tableRow( t[k+1] , ut , xt ,vt , ct, electricityPrice(t[k+1]));
            else
                output << t[k+1] <<" "<< ut <<" "<< xt << " "<< vt << " "<< ct << " " << electricityPrice(t[k+1])<< endl;
        }
        return 0;
    }
    
}

