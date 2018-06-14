#pragma once

namespace MSC_PROJECTS
{
// n is number of steps
// t_0 = a
// t_n = b
// t_i = a + ih where h=(b-a)/n 
// x_0 = alpha
// this will solve to find x_n
template<class F>
double eulersMethodTemplate(int n,double a,double b,
			    double alpha,const F &f)
{
  // local variables
  double h,t,x;
  // intialise values
  h=(b-a)/(double)(n);
  t=a;
  x=alpha;
  // implement Euler's method
  for(int i=0;i<n;i++)
  {
    t = a + i*h; // update value of t to t_i
    x = x + h*f(x,t); // update x to x_{i+1}
  }
  return x; // returns x_n
}

// n is number of steps
// t_0 = a
// t_n = b
// t_i = a + ih where h=(b-a)/n 
// x_0 = alpha
// this will solve to find x_n
template<class SCALAR,class VECTOR,class FUNCTION>
VECTOR midpointMethodTemplate(int n,SCALAR a,SCALAR b,
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
    VECTOR xMid = x + 0.5*h*f(x,t); // update x to x_{i+1}
    x = x + h*f(xMid,t+0.5*h); 
  }
  return x; // returns x_n
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

}