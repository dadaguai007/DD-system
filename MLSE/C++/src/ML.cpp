#include "ML.hpp"

ML::ML(uint8_t len_, int8_t CIR_len_)
{
    if (len_ > 30)
        throw runtime_error("possible input len " + to_string(len_) + " is incorrect. it should be less then 31");
    len                  = len_;
    CIR_len              = CIR_len_;
    conv_len             = len + CIR_len - 1;
    num_states           = 1 << len;
    states_complex_float = new complex<float>[num_states * len];
    states               = new int8_t[num_states * len];
    convolution_res      = new complex<float>[len + CIR_len - 1];
    gen_states();
}

ML::~ML()
{
    delete[] states;
    delete[] states_complex_float;
    delete[] convolution_res;
}

void ML::gen_states()
{
    for (size_t i = 0; i < num_states; i++) {
        int idx = i * len;
        for (int8_t j = len - 1; j >= 0; j--) {
            int8_t cur_bit                          = (i >> j) & 1;
            states[idx + len - j - 1]               = cur_bit == 0 ? -1 : 1;
            states_complex_float[idx + len - j - 1] = {static_cast<float>(states[idx + len - j - 1]), 0};
        }
    }
}

void ML::run(complex<float> *sig, complex<float> *CIR, int8_t *equilised_bits)
{
    if (sig == nullptr)
        throw runtime_error("sig is wrong");
    if (CIR == nullptr)
        throw runtime_error("CIR is wrong");
    complex<float> CIR_inverted[CIR_len];
    for (int i = 0; i < CIR_len; i++)
        CIR_inverted[i] = CIR[CIR_len - i - 1];
    float min_err = FLT_MAX;
    int   min_idx = 0;
    int   count   = 0;
    for (int i = 0; i < num_states; i++) {
        conv(&states_complex_float[i * len], CIR_inverted);
        float cur_err = diff(sig);
        if (cur_err < min_err) {
            min_idx = i;
            min_err = cur_err;
            // cout << cur_err << endl;
            count++;
        }
    }
    // if (count == 0)
    //     cout << "asnaljna " << count << endl;
    memcpy(equilised_bits, &states[min_idx * len], len * sizeof(states[0]));
}

float ML::diff(complex<float> *sig)
{
    float res = 0;
    for (int i = 0; i < conv_len; i++)
        res += abs(convolution_res[i] - sig[i]);
    return res;
}

void ML::conv(complex<float> *sig, complex<float> *CIR_inverted)
{
    if (CIR_inverted == nullptr)
        throw runtime_error("CIR_inverted is wrong");
    for (int i = 0; i < conv_len; i++) {
        if (i < CIR_len - 1)
            convolution_res[i] = inner_prod(sig, &CIR_inverted[CIR_len - i - 1], i + 1);
        else if (i >= len)
            convolution_res[i] = inner_prod(&sig[i - CIR_len + 1], CIR_inverted, CIR_len - (i - len + 1));
        else
            convolution_res[i] = inner_prod(&sig[i - CIR_len + 1], CIR_inverted, CIR_len);
    }
}

complex<float> ML::inner_prod(complex<float> *vec1, complex<float> *vec2, int vec_len)
{
    complex<float> res = 0;
    for (int i = 0; i < vec_len; i++)
        res += vec1[i] * vec2[i];
    return res;
}

complex<float> *ML::get_states()
{
    return states_complex_float;
}