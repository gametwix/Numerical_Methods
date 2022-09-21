#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include "../../Lab1/LB1.1/square_matrix.h"

std::vector<double> get_a(const std::vector<double> &Xi, const std::vector<double> &Fi, int n){
    SquareMatrix to_solve(n + 1);
    Vector<double> vec(n + 1);
    std::vector<double> a(n + 1);
    for(int k = 0; k < n + 1; ++k){
        for(int i = 0; i < n + 1; ++i){
            double sum = 0;
            for(int j = 0; j < Xi.size(); ++j){
                sum += std::pow(Xi[j],i+k);
            }
            to_solve(k,i) = sum;
        }
        double sum = 0;
        for(int j = 0; j < Xi.size(); ++j){
            sum += Fi[j]*std::pow(Xi[j],k);
        }
        vec(k) = sum; 
    }
    to_solve.find_LU();
    Vector<double> ans = to_solve.find_root(vec);
    for(int i = 0; i < ans.get_size(); ++i){
        a[i] = ans(i);
    }
    return a;
}

double F(double x, const std::vector<double> &a){
    double ans = 0;
    for(int i = 0; i < a.size(); ++i){
        ans += a[i]*std::pow(x,i);
    }
    return ans;
}

double find_eps(std::vector<double> Pred, std::vector<double> Fi){
    double ans = 0;
    for(int i = 0; i < Pred.size(); ++i){
        ans += std::pow(Pred[i] - Fi[i],2);
    }
    return ans;
}

int main(){
    //std::vector<double> Xi = {0.0, 1.7, 3.4, 5.1, 6.8, 8.5};
    //std::vector<double> Fi = {0.0, 1.3038, 1.8439, 2.2583, 2.6077, 2.9155};
    std::vector<double> Xi = {-0.9, 0.0, 0.9, 1.8, 2.7, 3.6};
    std::vector<double> Fi = {-1.2689, 0.0, 1.2689, 2.6541, 4.4856, 9.9138};
    std::vector<double> a1 = get_a(Xi,Fi,1);
    std::vector<double> a2 = get_a(Xi,Fi,2);

    std::cout << std::fixed << std::setprecision(4);

    //Print input data
    std::cout << " i ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << "   " << i << "   ";
    }
    std::cout << std::endl;
    std::cout << " xᵢ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << std::setw(7) << Xi[i];
    }
    std::cout << std::endl;
    std::cout << " yᵢ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << std::setw(7) << Fi[i];
    }
    std::cout << std::endl << std::endl;


    //Count F1(x) and F2(x)
    std::vector<double> F1(Xi.size());
    std::vector<double> F2(Xi.size());

    for(int i = 0; i < Xi.size(); ++i){
        F1[i] = F(Xi[i],a1);
    }
    for(int i = 0; i < Xi.size(); ++i){
        F2[i] = F(Xi[i],a2);
    }


    //Print F1
    std::cout << " i ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << "   " << i << "   ";
    }
    std::cout << std::endl;
    std::cout << " xᵢ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << std::setw(7) << Xi[i];
    }
    std::cout << std::endl;
    std::cout << " F₁";
    for(int i = 0; i < F1.size(); ++i){
        std::cout << std::setw(7) << F1[i];
    }
    std::cout << std::endl;

    std::cout << "Eps_F1 = " << find_eps(F1,Fi) << std::endl << std::endl;


    std::cout << " i ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << "   " << i << "   ";
    }
    std::cout << std::endl;
    std::cout << " xᵢ";
    for(size_t i = 0; i < 6; ++i){
        std::cout << std::setw(7) << Xi[i];
    }
    std::cout << std::endl;
    std::cout << " F₂";
    for(int i = 0; i < F2.size(); ++i){
        std::cout << std::setw(7) << F2[i];
    }
    std::cout << std::endl;

    std::cout << "Eps_F2 = " << find_eps(F2,Fi) << std::endl;

    std::ofstream out_file;
    out_file.open("points.txt");
    if(out_file.is_open()){
        for(int i = 0; i < Xi.size(); ++i){
            out_file << Xi[i] << " ";
        }
        out_file << std::endl;
        for(int i = 0; i < Fi.size(); ++i){
            out_file << Fi[i] << " ";
        }
        out_file << std::endl;
        for(int i = 0; i < F1.size(); ++i){
            out_file << F1[i] << " ";
        }
        out_file << std::endl;
        for(int i = 0; i < F2.size(); ++i){
            out_file << F2[i] << " ";
        }
        out_file << std::endl;
        out_file.close();
    }else{
        std::cout << "Error open output file" << std::endl;
    }

    return 0;
}