#include "mex.h"

void mexFunction(int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[])
{
    double* smallbins = mxGetPr(prhs[0]);
    double* P = mxGetPr(prhs[1]);
    unsigned int* sbi = mxGetPr(prhs[2]);
    unsigned int k = mxGetScalar(prhs[3]);
    int res = mxGetM(prhs[1]);
    unsigned int* sbiend = sbi + res;
    for ( ; sbi < sbiend; sbi++ ) {
        double* cur = smallbins + *sbi - k;
        double* endline = cur + res*k;
        for ( ; cur < endline; *P++ = *(cur+=k) );
    }
}