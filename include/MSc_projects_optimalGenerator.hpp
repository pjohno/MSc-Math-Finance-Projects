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
        typedef MATH60082::MVector Vector;
        typedef MATH60082::MMatrix Matrix;
        
        // grid for the state variable
        Vector x;
        
        // value function on grid points
        Vector v;
        // value function on grid points
        Matrix control;
        
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
        
        typedef MATH60082::MVector Vector;
        typedef MATH60082::MMatrix Matrix;
        typedef std::vector<int> CV;
        
        std::function<double(double)> EP;
        
        double tol;
        
        Generator(int n,int m,double T,double tol);
        
        std::vector<GeneratorState> G;
        Vector t;
        
        int solve();
        
        int outputPath(double x0,int u0);
        
        template <class T>
        int setupState(int state,int n,int m,double tol,const T &params);
        
    };
    
}

