#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>
#include "math60082_mVector.hpp"
#include "math60082_mMatrix.hpp"
#include "math60082_lagrange_interp.hpp"

namespace MSC_PROJECTS
{
    
    
    struct GeneratorState
    {
        // grid for the state variable
        MATH60082::MVector x;
        
        // value function on grid points
        MATH60082::MVector v;
        // value function on grid points
        MATH60082::MMatrix control;
        
        // value function operator
        double operator()(double x0){return MATH60082::lagrangeInterpolation(v.returnArray(),x.returnArray(),x0,v.size(),2);}
        
        // dynamics of the state
        double f(double,double) const;
        std::function<double(double)> cMin;
        std::function<double(double)> cMax;
        std::function<int(double)> eta;
        
        // cost function
        std::function<double(double,double,double)> Gamma;
        
    };
    
    struct GeneratorDefaultOff{
        double xMin;
        double xMax;
        double cMin;
        double cMax;
        double cost;
    };
    
    struct GeneratorDefaultSyncing{
        double xMin;
        double xMax;
        double cMin;
        double cMax;
        double cost;
    };
    
    struct GeneratorDefaultDeSyncing{
        double xMin;
        double xMax;
        double cMin;
        double cMax;
        double cost;
    };
    
    struct GeneratorDefaultOn{
        double xMin;
        double xMax;
        double cMin;
        double cMax;
        std::function<double(double,double,double)> cost;
    };
    
    struct Generator
    {
        
        std::function<double(double)> EP;
        
        double tol;
        
        int defaultSetup(int n,int m,double T,double tol);
        
        std::vector<GeneratorState> G;
        MATH60082::MVector t;
        
        int solve();
        
        int outputPath(double x0,int u0,std::ostream& output=std::cout,bool toMarkup=true);
        
        template <class T>
        int setupState(int state,int n,int m,double tol,const T &params);
        
        GeneratorState& operator[](unsigned int U){return G[U];}
        
    };
    
}

