#include "mex.h"
// #include "class_handle.hpp"
#include "C++/src/ML.hpp"
#include <string>


static ML* obj = nullptr;
static int8_t bit_num;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    

    if (!mxIsChar(prhs[0])) {
        mexErrMsgTxt("first input must be a char array");
    }
    char* command = mxArrayToString(prhs[0]);

    if (strcmp(command, "new") == 0) {
        if (obj != nullptr) {
            delete obj;
        }
        bit_num = mxGetScalar(prhs[1]);
        int8_t CIR_len = mxGetScalar(prhs[2]);
        obj = new ML(bit_num, CIR_len);
    }else if (strcmp(command, "run") == 0) {
        if (obj != nullptr) {
            std::complex<float>* sig = (std::complex<float>*)mxGetData(prhs[1]);
            std::complex<float>* CIR = (std::complex<float>*)mxGetData(prhs[2]);
            std::int8_t equilised_bits[bit_num];
            obj->run(sig, CIR, equilised_bits);
            plhs[0] = mxCreateNumericMatrix(1, bit_num, mxINT8_CLASS, mxREAL);
            std::memcpy(mxGetData(plhs[0]), equilised_bits, bit_num * sizeof(int8_t));
        } else {
            mexErrMsgTxt("Object not created.");
        }
    } else {
        mexErrMsgTxt("Unknown command.");
    }

    mxFree(command);
}
