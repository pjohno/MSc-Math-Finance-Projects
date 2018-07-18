

namespace MSC_PROJECTS
{
    /** @brief Newton secant method 
     * @details Find the root of a function using a newton secant method.
     * @param xold initial guess \f$ x_0 \f$
     * @param xcurr input initial guess \f$ x_1 \f$, on return xcurr is the root
     * @param maxiter is the max no of iterations
     * @param tol is the accuracy of the root
     * @return 0 will be returned if no root is found, 1 if it is successful
     */
    template<class T>
    static int findroot (double xold , double& xcurr , int maxiter , double tol , const T& f){ 
        // find root
        int iter ;
        // x_{n} is xcurr and x_{n-1} is xold
        // temporary variable to store x_{n+1}
        double xnew,fold=f(xold),fcurr;
        // algorithm given by
        //   x_{n+1} = x_n - f(x_n)(x_{n}-x_{n-1}) / (f(x_n)-f(x_{n-1}))
        for ( iter =0; iter <maxiter ; iter++){
            fcurr=f(xcurr);
            xnew = xcurr - fcurr * ( xcurr - xold )/ ( fcurr - fold ) ;
            // std::cout << xnew << " " << f(xcurr) << "\n";
            fold = fcurr;
            xold = xcurr;
            xcurr = xnew;
            if ( f(xcurr) < tol && f(xcurr) > -tol )break ;
        }
        if ( iter==maxiter ) return 0 ;
        else return 1 ;
    }
    
}
