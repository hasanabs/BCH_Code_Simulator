// Written by Hasan Albinsaid <hasanalbinsaid@hotmail.com>
//For more information: https://github.com/hasanabs
#include "fungsi.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <vector>
#include <math.h>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);
normal_distribution<double> distribution (0.0,1.0);

vector<int> ECC::mul_poly_bin(vector<int> &x, vector<int> &y){
    vector<int> result;
    for (int i=0;i<(int)x.size()+(int)y.size()-1;i++) result.push_back(0);
    for (int i=0;i<(int)x.size();i++){
        rotate(result.begin(),result.begin()+result.size()-1,result.end()); //cyclic rotation. Minus (-) denote delay or shift to right
        for (int j=0;j<(int)y.size();j++){
            result[j]=(result[j]+y[j]*x[(int)x.size()-i-1])%2;
        }
    }
    while(result[(int)result.size()-1]==0)result.pop_back();
    return result;
}

vector<int> ECC::cyclic_mul_poly_bin(vector<vector<int>> &a_table, vector<int> &x, vector<int> &y){
    vector<int> poly_result = ECC::mul_poly_bin(x,y);
    vector<int> result;
    result.insert(result.end(),a_table.at(0).begin(),a_table.at(0).end());
    for (int i=0; i<(int)poly_result.size(); i++){
        if (poly_result[i]!=0){
            vector<int> temp_result = xor_poly_bin(a_table.at(i+1), result);
            result.swap(temp_result);
        }
    }
    return result;
}

vector<int> ECC::cyclic_div_poly_bin(vector<vector<int>> &a_table, vector<int> &x, vector<int> &y){
    int diferent = ECC::index_table_finder(a_table,x)-ECC::index_table_finder(a_table,y);
    if (diferent>=0) {return a_table.at(diferent+1);}
    else {return a_table.at((int)a_table.size()+diferent);}
}

vector<int> ECC::mod_poly_bin(vector<int> &x, vector<int> &y){
    vector<int> result;
    result.insert(result.end(),x.begin(),x.end()); //Put "x vector" to end of "result vector"
    vector<int> deduction;
    vector<int> temp;

    while (result[(int)result.size()-1]==0)   result.pop_back();
    while ((int)result.size()>=(int)y.size()){
        deduction.clear();
        deduction.insert(deduction.end(),y.begin(),y.end());
        while((int)result.size()>(int)deduction.size()) deduction.insert(deduction.begin(), 0); //push_front
        temp.clear();   temp = ECC::xor_poly_bin(result, deduction);
        result.clear(); for (int i=0;i<(int)temp.size();i++) result.push_back(temp[i]);
        while (result[(int)result.size()-1]==0)   result.pop_back();
    }
    return result;
}

vector<int> ECC::div_poly_bin(vector<int> &x, vector<int> &y){
    vector<int> result;
    vector<int> devidend=x; //Put "x vector" to end of "devidend vector"

    while ((int)devidend.size()>=(int)y.size()){
        if (devidend[(int)devidend.size()-1]==0){
            result.insert(result.begin(),0);
            devidend.pop_back();
        }
        else{
            result.insert(result.begin(),1);
            vector<int> divisor=y;
            while((int)devidend.size()>(int)divisor.size()) divisor.insert(divisor.begin(), 0); //push_front
            devidend = ECC::xor_poly_bin(devidend, divisor);
            devidend.pop_back();
        }
    }
    return result;
}

vector<int> ECC::xor_poly_bin(vector<int> &a, vector<int> &b){
    vector<int> result;
    for (int i=0;i<(int)a.size();i++) result.push_back((a[i]+b[i])%2);
    return result;
}

vector<vector<int>> ECC::primitive_poly_table_bin(int GF,vector<int> &primitive_poly){
    vector<vector<int>> matrik;
    vector<int> vect;
    for(int i=0;i<GF;i++) vect.push_back(0); // add 0,0,..,0
    matrik.push_back(vect);
    vect.clear();    vect.push_back(1);    for(int i=1;i<GF;i++) vect.push_back(0); // add 1,0,..,0
    matrik.push_back(vect);
    for (int i=2;i<(int)pow(2,GF);i++){
        vect.clear();
        if (matrik[i-1][GF-1]==1){
            for (int j=0;j<GF;j++){
                if (j==0) vect.push_back((primitive_poly[j]+0)%2);
                if (j>0) vect.push_back((primitive_poly[j]+matrik[i-1][j-1])%2);
            }
        }
        else{
            for (int j=0;j<GF;j++){
                if (j==0) vect.push_back(0);
                if (j>0) vect.push_back(matrik[i-1][j-1]);
            }
        }
        matrik.push_back(vect);
    }
    return matrik;
}

vector<vector<vector<int>>> ECC::BCH_parity(int t,int GF, vector<vector<int>> &a_table){
    vector<vector<vector<int>>> result;
    vector<vector<int>> temp1;
    vector<int> temp2, temp2_static;
    for (int i=0; i<2*t; i++){
        temp1.clear();
        for (int j=1; j<(int)pow(2,GF); j++){
            temp2.clear(); temp2_static.clear();
            temp2 = a_table.at(j); temp2_static = a_table.at(j);
            for (int m=0; m<i; m++){
                temp2 = ECC::cyclic_mul_poly_bin(a_table,temp2, temp2_static);
            }
            temp1.push_back(temp2);
        }
        result.push_back(temp1);
    }
    return result;
}

vector<vector<int>> ECC::syndrome_BCH(vector<int> &receive, vector<vector<int>> &a_table, vector<vector<vector<int>>> &H){
    vector<vector<int>> result;
    vector<int> temp_res, temp_plus;
    for (int i=0; i<(int)H.size(); i++){
    temp_plus.clear();
    temp_plus.insert(temp_plus.end(),a_table.at(0).begin(),a_table.at(0).end());
        for (int j=0; j<(int)receive.size(); j++){
            if (receive[j]!=0){
                temp_plus=ECC::xor_poly_bin(temp_plus,H.at(i).at(j));
            }
        }
       result.push_back(temp_plus);
    }
    return result;
}

vector<int> ECC::encode_BCH(vector<int> &data, vector<int> &gx, int n_k){
    vector<int> result=data;
    for(int i=0; i<n_k; i++) result.insert(result.begin(), {0}); //push_front
    vector<int> reminder = ECC::mod_poly_bin(result,gx);reminder.resize(n_k);
    for(int i=0; i<n_k; i++) result.at(i)=reminder.at(i); //shift left
    return result;
}

vector<vector<int>> ECC::decode_BCH(vector<vector<int>> &a_table, vector<vector<int>> syndrome){
    vector<int> zeros;
    zeros.insert(zeros.end(),a_table.at(0).begin(),a_table.at(0).end());
    /*-- Initialization table --*/
    int rho=0; int trace=0;
    vector<int> poly_init {1};
    for (int i=1;i<(int)syndrome.at(0).size();i++) poly_init.push_back(0);
    vector<vector<vector<int>>> sigma {{poly_init}};  sigma.push_back({{poly_init}});
    vector<vector<int>> du {poly_init}; du.insert(du.end(),syndrome.at(0));
    vector<int> lu {0,0};
    vector<int> u_lu {-1,0};
    /*-- Initialization table --*/

    for (int i=2; i<(int)syndrome.size()+2; i++){
        if (accumulate(du.at(i-1).begin(),du.at(i-1).end(),0)==0){ //if du=0 next sigma will equal with before
            sigma.push_back({{sigma.at(i-1)}}); //use old value of sigma
        }
        else{ //if du!=0
            sigma.push_back({{sigma.at(i-1)}});
            for (trace=i-2; trace>=0; trace--) if (accumulate(du.at(trace).begin(),du.at(trace).end(),0)!=0) break;
            rho = trace; //u[trace];
            vector<int> du_dp = ECC::cyclic_div_poly_bin(a_table,du.at(i-1),du.at(rho));
            vector<vector<int>> x_times_sigmarho = sigma.at(rho); for(int j=0; j<(i-1)-rho; j++) x_times_sigmarho.insert(x_times_sigmarho.begin(), zeros); //push_front

            for (int k=0; k<(int)x_times_sigmarho.size();k++){ //du/dp*x*sigma
                if (accumulate(x_times_sigmarho.at(k).begin(), x_times_sigmarho.at(k).end(), 0)!=0){
                    x_times_sigmarho.at(k)=ECC::cyclic_mul_poly_bin(a_table, x_times_sigmarho.at(k), du_dp);
                }
            }
            if ((int)x_times_sigmarho.size()>(int)sigma.at(i).size()){
                for(int l=0;l<(int)x_times_sigmarho.size()-(int)sigma.at(i).size();l++) sigma.at(i).insert(sigma.at(i).end(), zeros);
            }
            else if ((int)x_times_sigmarho.size()<(int)sigma.at(i).size()){
                for(int l=0;l<(int)sigma.at(i).size()-(int)x_times_sigmarho.size();l++) x_times_sigmarho.insert(x_times_sigmarho.end(), zeros);
            }
            for (int k=0; k<(int)sigma.at(i).size();k++){
                sigma.at(i).at(k)= ECC::xor_poly_bin(sigma.at(i).at(k),x_times_sigmarho.at(k)); //use old value of sigma + corrction termk
            }
        }
        if (i==(int)syndrome.size()+1) break;
        lu.push_back((int)sigma.at(i).size()-1);
        u_lu.push_back(i-1-lu.at(i));
        //Update new du
        du.push_back(syndrome.at(i-1)); //s_(u+1)
        for (int k=0; k<lu[i];k++){
            vector<int> du_multi = ECC::cyclic_mul_poly_bin(a_table, sigma.at(i).at(k+1),syndrome.at(i-k-2));
            du.at(i) = ECC::xor_poly_bin(du.at(i),du_multi);
        }
    }
    return sigma.at((int)syndrome.size()+1);
}

vector<int> ECC::error_loc_BCH(vector<vector<int>> &a_table, vector<vector<int>> &sigma_x){
    vector<int> result;
    for (int i=1;i<(int)a_table.size();i++){
        vector<int> temp = {1};
        temp.insert(temp.end(),a_table.at(0).begin(),a_table.at(0).end()-1);//push 1,0,...,0
        for (int j=1;j<(int)sigma_x.size();j++){
            vector<int> a_static=a_table.at(i);
            for(int k=1; k<j; k++){ //power of x (x^j)
                a_static=ECC::cyclic_mul_poly_bin(a_table, a_static,a_table.at(i));
            }
            a_static=ECC::cyclic_mul_poly_bin(a_table, a_static,sigma_x.at(j));
            temp=ECC::xor_poly_bin(temp,a_static);
        }
        if (accumulate(temp.begin(),temp.end(),0)==0) result.push_back(1);
        if (accumulate(temp.begin(),temp.end(),0)!=0) result.push_back(0);
    }
    reverse(result.begin(),result.end()); //Beta => inverse of Alfa
    rotate(result.begin(),result.begin()+result.size()-1,result.end()); //Beta => inverse of Alfa
    return  result;
}

int ECC::index_table_finder(vector<vector<int>> &a_table, vector<int> &x){
    int i;
    for (i=0; i<(int)a_table.size(); i++){
        if (a_table.at(i)==x) break;
    }
    return i;
}

void ECC::print_GF_table(vector<vector<int>> &a_table){
    cout << "     Poly Table    "<<endl;
    cout <<"GF(2^"<<(int)a_table.at(0).size()<< ") |<LSB- -MSB>" << endl;
    for (int i=0;i<(int)a_table.size();i++){ //Print the binary representation GF(2^x)
        if (i==0) cout << i << "\t| ";
        if (i>0) cout << "a^" << i-1 << "\t| ";
        for (int j=0;j<(int)a_table.at(0).size();j++) cout << a_table[i][j] << " ";
        cout << endl;
    }
}

void ECC::init_seed(){
    srand(time(0));
}

vector<int> ECC::BPSK_modulation(vector<int> &x){
    vector <int> result;
    for (int i=0; i<(int)x.size();i++){
        result.push_back(x.at(i)*2-1);
    }
    return  result;
}

vector<int> ECC::BPSK_demodulation(vector<float> &x){
    vector <int> result;
    for (int i=0; i<(int)x.size();i++){
        x.at(i)>0 ? result.push_back(1) : result.push_back(0);
    }
    return  result;
}

vector<float> ECC::add_AWGN_noise(vector<int> &x, float R, int snr){
    vector <float> result;
    for (int i=0; i<(int)x.size();i++){
        float sigma_square = 1.0/(2*R*pow(10,(float)snr/10));
        float noise = sqrt(sigma_square) * distribution(generator);
        result.push_back(x.at(i)+noise);
    }
    return  result;
}

vector<int> ECC::random_bin_data(int length){
    vector<int> result;
    for (int i=0; i<length; i++) result.push_back(rand()%2);
    return  result;
}

vector<int> ECC::repair_bit(vector<int> &x, vector<int> &y){
    vector <int> result;
    for (int i=0; i<(int)x.size();i++){
        result.push_back((x.at(i)+y.at(i))%2);
    }
    return  result;
}
