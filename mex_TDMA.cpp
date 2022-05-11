/**************************************************************************
*	TDMA - solve tridiagonal linear system using Tridiagonal matrix 
*       algorithm. (no LR decomposition!)
*
*	Usage: u = mex_TDMA(a, b, c, f)
*
*	Input:
*       'a', 'b', 'c' - main, upper and lower diagonals of a square  
*           tridiagonal matrix A of size NxN.
*       'f' - vector of size Nx1 or 1xN.
*	Output:
*       'u' - solution of the system Au = f.
*
*   Compile using 'mex -largeArrayDims mex_TDMA.cpp' in Matlab command line.
*
*	Copyright (c) 2022 Alexey Tarasov, based on 
*       https://www.mathworks.com/matlabcentral/fileexchange/38640-implementation-of-thomas-algorithm-mex
**************************************************************************/
#include "mex.h"
#include <vector>
#include <math.h>

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *a, *b, *c, *f;
    std::vector<double> alpha, beta;
    double *u;
    size_t N;
    
    /* Check for proper number of input and output arguments */
    if (nrhs != 4) {
        mexErrMsgTxt("Four input arguments required: a, b, c, f.");
    }
    if(nlhs > 4){
        mexErrMsgTxt("Too many output arguments.");
    }

    N = mxGetNumberOfElements(prhs[0]);
     
    if ( (mxGetNumberOfElements(prhs[1]) != N - 1) || 
            (mxGetNumberOfElements(prhs[2]) != N - 1) ||
            (mxGetNumberOfElements(prhs[3]) != N) ) {
        mexErrMsgTxt("The following should hold: |b| = |c| = |a| - 1, |f| = |a|.");
    }
    
    a = (double *)mxGetData(prhs[0]);
    b = (double *)mxGetData(prhs[1]);
    c = (double *)mxGetData(prhs[2]);
    f = (double *)mxGetData(prhs[3]);
    
    alpha.resize(N, 0);
    beta.resize(N, 0);
    
    alpha[0] = -b[0]/a[0];
    beta[0] = f[0]/a[0];
    for (size_t i = 1; i < N-1; i++) {
        alpha[i] = -b[i] / ( c[i-1]*alpha[i-1] + a[i] );
        beta[i] = (f[i] - c[i-1]*beta[i-1])/(c[i-1]*alpha[i-1] + a[i]); 
    }
    beta[N-1] = (f[N-1] - c[N-2]*beta[N-2])/( c[N-2]*alpha[N-2] + a[N-1] );
    //beta[N-1] = (f[N-1] - c[N-2]*beta[N-2])/( c[N-2]*alpha(N-1) + A(N,N) );
    
    plhs[0] = mxCreateDoubleMatrix((mwSize)N, 1, mxREAL);
    u = (double *)mxGetData(plhs[0]);
    u[N-1] = beta[N-1];
    //u[N-1] = beta[N-1];
    for (size_t i = N-2; i > 0; i--) {
        u[i] = alpha[i]*u[i+1] + beta[i];
        //u[i] = beta[i];
    }
    u[0] = alpha[0]*u[1] + beta[0];   
    //u[0] = beta[0];
}
