#include "mex.h"
#include <cmath>
#include <cstring>
#include <fstream>

void mexFunction(int nlhs, mxArray *plhs[ ], int nrhs, const mxArray *prhs[ ]) {
    // The order of input
    // (1) p matrix, (2) NN matrix
    // Src samples: n x d array
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
            output[j*nSrc + i] =  ptr1[j*nSrc + i] + ptr2[j*nSrc + i];
        }
    }
    
    return;
}