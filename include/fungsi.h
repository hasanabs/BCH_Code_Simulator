// Written by Hasan Albinsaid <hasanalbinsaid@hotmail.com>
//For more information: https://github.com/hasanabs
//You can read all the description of each function in this file
#ifndef FUNGSI_H
#define FUNGSI_H
#include <vector>
using namespace std;


class ECC{
    public:
        //To multiply polynomial in GF(2), Note: this is not cyclic. [input: 2 polinomial represent in binary, MSB is in right side]
        vector<int> mul_poly_bin(vector<int> &x, vector<int> &y);

        //To multiply polynomial in cyclic GF(2) according to GF table. [input: GF table; 2 polinomial represent in binary, MSB is in right side]
        vector<int> cyclic_mul_poly_bin(vector<vector<int>> &a_table, vector<int> &x, vector<int> &y);

        //To devide polynomial in cyclic GF(2) according to GF table. [input: GF table; 2 polinomial represent in binary, MSB is in right side]
        vector<int> cyclic_div_poly_bin(vector<vector<int>> &a_table, vector<int> &x, vector<int> &y);

        //To devide polynomial in GF(2) according to GF table. [input: 2 polinomial represent in binary, MSB is in right side]
        vector<int> mod_poly_bin(vector<int> &x, vector<int> &y);

        //To devide polynomial in GF(2) according to GF table. [input: 2 polinomial represent in binary, MSB is in right side]
        vector<int> div_poly_bin(vector<int> &x, vector<int> &y);

        //To do adition or substraction in GF(2), Note: lenght of bit should equal. [input: 2 polinomial represent in binary, MSB is in right side]
        vector<int> xor_poly_bin(vector<int> &x, vector<int> &y);

        //To generate GF table according to input GF and primitive polynomial. [input: m degree of GF(2^m); primitive polinomial represent in binary, MSB is in right side]
        vector<vector<int>> primitive_poly_table_bin(int GF,vector<int> &primitive_poly);

        //To generate parity check matrix for BCH code. [input: t error-correction; m degree of GF(2^m); GF table]
        vector<vector<vector<int>>> BCH_parity(int t,int GF, vector<vector<int>> &a_table);

        //To get syndrom of BCH code. [input: receive polynomial represent in binary, MSB is in right side; GF table; BCH parity check matrix]
        vector<vector<int>> syndrome_BCH(vector<int> &receive, vector<vector<int>> &a_table, vector<vector<vector<int>>> &H);

        //To encode systematic BCH code. [input: GF table; syndrom of BCH code]
        vector<int> encode_BCH(vector<int> &data, vector<int> &gx, int n_k);

        //To get the value of sigma(x) in BCH code. [input: GF table; syndrom of BCH code]
        vector<vector<int>> decode_BCH(vector<vector<int>> &a_table, vector<vector<int>> syndrome);

        //To get the error bit location in receiver BCH code. [input: GF table; sigma(c) of BCH code]
        vector<int> error_loc_BCH(vector<vector<int>> &a_table, vector<vector<int>> &sigma_x);

        //Another helper function
        int index_table_finder(vector<vector<int>> &a_table, vector<int> &x);
        void print_GF_table(vector<vector<int>> &a_table);
        void init_seed();
        vector<float> add_AWGN_noise(vector<int> &x, float R, int snr);
        vector<int> random_bin_data(int length);
        vector<int> BPSK_modulation(vector<int> &x);
        vector<int> BPSK_demodulation(vector<float> &x);
        vector<int> repair_bit(vector<int> &x, vector<int> &y);
};

#endif // FUNGSI_H
