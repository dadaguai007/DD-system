#include <math.h>
#include <iostream>
#include <complex>
#include <float.h>
#include <cstring>
using namespace std;

class ML
{
public:
    ML(uint8_t len_, int8_t CIR_len_);
    ~ML();
    complex<float> *get_states();
    void            run(complex<float> *sig, complex<float> *CIR, int8_t *equilised_bits);

private:
    uint8_t         len;
    int8_t          CIR_len;
    int8_t          conv_len;
    int             num_states;
    complex<float> *convolution_res;
    complex<float> *states_complex_float;
    int8_t         *states;
    void            gen_states();
    void            conv(complex<float> *sig, complex<float> *CIR);
    complex<float>  inner_prod(complex<float> *vec1, complex<float> *vec2, int vec_len);
    float           diff(complex<float> *sig);
};