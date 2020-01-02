#include <iostream>
#include <numeric>
#include <math.h>
#include<vector>
#include <fungsi.h>

using namespace std;
ECC hasan;

int t=2;
int GF=4;
vector <int> primitive_poly = {1,1,0,0}; //x^4 = 1 + 1*x + 0*x^2 + 0*x^3
vector <int> g_x = {1,0,0,0,1,0,1,1,1}; //g(x)=1 + x^4 + x^6 + x^7 + x^8

int main()
{
    int snr_min=0, snr_max=11;
    long rep=10000000;
    hasan.init_seed(); //To make true random
    int n = pow(2,GF) -1;
    int k = pow(2,GF)-(int)g_x.size();
    float R = float(k)/(float)(pow(2,GF)-1);
    vector<float> error; long repeat;

    /*---------- Generate GF table and Parity check ----------*/
    vector<vector<int>> a_table = hasan.primitive_poly_table_bin(GF,primitive_poly); //Generate representation GF(2^x)
    vector<vector<vector<int>>> H = hasan.BCH_parity(t,GF,a_table); //Parity check matrix
    hasan.print_GF_table(a_table); cout <<endl<<endl;
    cout <<"|  SNR  ||\t   BER   \t|" <<endl;

    /*---------- main loop ----------*/
    for (int snr=snr_min; snr<=snr_max; snr++){
        long error_snr=0;
        for (repeat=0; repeat<rep; repeat++){
            vector <int> data = hasan.random_bin_data(k);
            vector <int> encode = hasan.encode_BCH(data,g_x,n-k);
            vector <int> transmit = hasan.BPSK_modulation(encode);
            /*---------- Transmit trough channel ----------*/
            vector <float> receive = hasan.add_AWGN_noise(transmit,R,snr);
            vector <int> receive_demod = hasan.BPSK_demodulation(receive);
            vector<vector<int>> syndrome =hasan.syndrome_BCH(receive_demod,a_table,H);
            vector<vector<int>> sigma_x = hasan.decode_BCH(a_table, syndrome);
            vector<int> error_loc = hasan.error_loc_BCH(a_table,sigma_x);
            vector<int> fix_error = hasan.repair_bit(receive_demod,error_loc); fix_error.shrink_to_fit();
            vector<int> decoded (fix_error.cbegin()+n-k,fix_error.cend());
            vector<int> different = hasan.xor_poly_bin(data,decoded);
            error_snr = error_snr + accumulate(different.begin(),different.end(),0);
        }
        error.push_back((float)error_snr/(float)k/(float)repeat);
        cout <<fixed; cout <<"   " <<snr;
        cout <<scientific; cout.precision(4); cout <<"\t\t " <<error.at(snr-snr_min) <<endl;
    }
}
