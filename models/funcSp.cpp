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

inline void funcSp(double mu, double nu, double * res){
    double z = 0;
    double temp = 0;
    res[0] = 0;
    res[1] = 0;
    for (int i = 0; i <= 7; i++){
        temp = pow(double(mu), double(i*nu))/pow(double(fakt(i)), double(nu));
        z += temp;
        res[0] += i*nu*temp/mu;
        res[1] += log(pow(double(mu), double(i))/fakt(i))*temp;
    }  
    
    res[0] /= z;
    res[1] /= z;
}

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {
    if(nrhs!=2) {
        mexErrMsgTxt("Two inputs are needed.\n");
    }
    int i = 0;
    int j = 0;
    double res[2];
    int nSrc = mxGetM(prhs[0]);
    int dSrc = mxGetN(prhs[0]);
    if( (nSrc != mxGetM(prhs[1])) || (dSrc != mxGetN(prhs[1]))){
        mexErrMsgTxt("Input matrices should have the same size.\n");
    }

    double *ptr1 = (double *)mxGetPr(prhs[0]);
    double *ptr2 = (double *)mxGetPr(prhs[1]);
    
    double *Gmu;
    plhs[0] = mxCreateDoubleMatrix(nSrc, dSrc, mxREAL);
    Gmu = mxGetPr(plhs[0]);
    double *Gnu;
    plhs[1] = mxCreateDoubleMatrix(nSrc, dSrc, mxREAL);
    Gnu = mxGetPr(plhs[1]);
    
    // Real program
    for (i=0; i<nSrc; i++){
        for (j=0; j<dSrc; j++){
            funcSp(ptr1[j*nSrc + i], ptr2[j*nSrc + i], res);
            Gmu[j*nSrc + i] =  res[0];
            Gnu[j*nSrc + i] =  res[1];
        }
    }
    
    return;
}