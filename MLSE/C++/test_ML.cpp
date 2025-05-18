#include <ML.hpp>
#include <iostream>
using namespace std;

int main()
{
    const int bit_len = 16;
    const int CIR_len = 2;
    const int N       = bit_len + CIR_len - 1;
    ML        equaliser(bit_len, CIR_len);
    float     rec_sig[N * 2] = {0.458922034450140, 0.536307469722595,  -0.376288483530962, 1.32092237374003,   -0.376288483530962, 1.32092237374003,
                                -1.29413255243124, 0.248307434294841,  1.29413255243124,   -0.248307434294841, -1.29413255243124,  0.248307434294841,
                                1.29413255243124,  -0.248307434294841, -1.29413255243124,  0.248307434294841,  1.29413255243124,   -0.248307434294841,
                                -1.29413255243124, 0.248307434294841,  1.29413255243124,   -0.248307434294841, -1.29413255243124,  0.248307434294841,
                                1.29413255243124,  -0.248307434294841, -0.376288483530962, 1.32092237374003,   -1.29413255243124,  0.248307434294841,
                                0.376288483530962, -1.32092237374003,  0.835210517981102,  -0.784614904017436};

    float  cir[CIR_len * 2] = {-0.458922034450140, -0.536307469722595, 0.835210517981102, -0.784614904017436};
    int8_t tx_bits[bit_len] = {0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1};
    int8_t rx_bits[bit_len];
    equaliser.run(reinterpret_cast<complex<float> *>(rec_sig), reinterpret_cast<complex<float> *>(cir), rx_bits);
    int8_t err = 0;
    for (int i = 0; i < bit_len; i++)
        err += abs(tx_bits[i] * 2 - 1 - rx_bits[i]);
    if (err == 0)
        cout << "TEST PASS" << endl;
    else {
        cout << "TEST FAIL! BER is " << static_cast<int>(err) / bit_len << endl;
    }
    return 0;
}