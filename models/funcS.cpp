#include "mex.h"
#include <cmath>
#include <cstring>

inline double fakt(int x){
    double out = 1;
    for (int i = 1; i <= x; i++){
        out *= i;
    }
    return out;
}

inline double funcS(double mu, double nu){
    double z = 0;
    for (int i = 0; i <= 7; i++){
        z += pow(double(mu), double(i*nu))/pow(double(fakt(i)), double(nu));
    }  
    
    return log(z);
}

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {
    if(nrhs!=2) {
        mexErrMsgTxt("Two inputs are needed.\n");
    }
    int i = 0;
    int j = 0;
    int nSrc = mxGetM(prhs[0]);
    int dSrc = mxGetN(prhs[0]);
    if( (nSrc != mxGetM(prhs[1])) || (dSrc != mxGetN(prhs[1]))){
        mexErrMsgTxt("Input matrices should have the same size.\n");
    }

    double *ptr1 = (double *)mxGetPr(prhs[0]);
    double *ptr2 = (double *)mxGetPr(prhs[1]);
    
    double *output;
    plhs[0] = mxCreateDoubleMatrix(nSrc, dSrc, mxREAL);
    output = mxGetPr(plhs[0]);
    
    // Real program
    for (i=0; i<nSrc; i++){
        for (j=0; j<dSrc; j++){
            output[j*nSrc + i] =  funcS(ptr1[j*nSrc + i], ptr2[j*nSrc + i]);
        }
    }
    
    return;
}