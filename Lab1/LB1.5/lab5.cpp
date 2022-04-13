#include <iostream>
#include <cmath>
#include <tuple>
#include <complex>
#include <string>
#include "../matrix.h"

std::string subs[10] = {"₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"};

template <class T>
T sign(T num){
    if(num > 0) return 1;
    else if(num == 0) return 0;
    else return -1;
}

template <class T>
Matrix<T> find_H(Matrix<T> mat, int step){
    Matrix<T> E(mat.get_row(),mat.get_row());
    Vector<T> V(mat.get_row());
    for(int i = 0; i < E.get_row(); ++i){
        E(i,i) = 1;

        if(i < step){
            V(i) = 0;
        }else if(i > step){
            V(i) = mat(i, step);
        }else{
            double sum = 0;
            for(int j = step; j < mat.get_row(); ++j){
                sum += mat(j, step)*mat(j, step);
            }
            sum = std::sqrt(sum);
            V(i) = mat(i,i) + sign(mat(i,i))*sum;
        }
    }

    Matrix<T> H = E - (V*V.transposition())*(2 / (V.transposition()*V)(0,0));
    return H;
}


template <class T>
std::pair<Matrix<T>,Matrix<T>> find_QR(Matrix<T> mat){
    Matrix<T> Q(mat.get_row(),mat.get_row());
    for(int i = 0; i < mat.get_row(); ++i){
        Q(i,i) = 1;
    }
    int steps = mat.get_row() - 1;
    for(int i = 0; i < steps; ++i){
        Matrix<T> H = find_H(mat,i);
        Q = Q*H;
        mat = H*mat;
    }
    return {Q,mat};
}

bool check_cols(Matrix<double> mat, double eps){
    //check eps
    bool check = true;
    for(int i = 0; i < mat.get_col(); ++i){
        //real number
        double sum = 0;
        for(int j = i + 1; j < mat.get_row();++j){
            sum += mat(j,i)*mat(j,i);
        }
        sum = sqrt(sum);
        if(sum <= eps){
            continue;
        }
        //complex number
        sum = 0;
        for(int j = i + 2; j < mat.get_row();++j){
            sum += mat(j,i)*mat(j,i);
        }
        sum = sqrt(sum);
        if(sum <= eps){
            continue;
        }
        check = false;
    }
    return check;
}

bool IsRealEigenv(Matrix<double> mat, int col, double eps){
    double sum = 0;
    for(int j = col + 1; j < mat.get_row();++j){
        sum += mat(j,col)*mat(j,col);
    }
    sum = sqrt(sum);
    return sum <= eps;
}

std::pair<std::complex<double>,std::complex<double>> find_pair_compl_num(Matrix<double> mat,int col){
    std::complex<double> in1,in2,in3;
    in1 = mat(col,col);
    in2 = mat(col + 1,col + 1);
    in3 = mat(col + 1,col) * mat(col,col + 1);
    std::complex<double> b,c,d;
    b = in1+in2;
    c = in1*in2 - in3;
    d = b*b - 4.*c;
    std::complex<double> ans1,ans2;
    ans1 = (b + sqrt(d))/2.;
    ans2 = (b - sqrt(d))/2.;
    return {ans1,ans2};
}

std::vector<std::complex<double>> find_eigenvalues(Matrix<double> mat, double eps){
    Matrix<double> last = mat;
    std::vector<std::complex<double>> answer(mat.get_col());
    while(true){
        Matrix<double> Q,R;
        std::tie(Q,R) = find_QR(mat);
        last = mat;
        mat = R*Q;
        if(check_cols(mat,eps)){
            bool check = true;
            for(int i = 0; i < mat.get_col();++i){
                if(IsRealEigenv(mat,i,eps)){
                    answer[i] = mat(i,i);
                }else{
                    if(IsRealEigenv(mat,i + 1,eps)){
                        std::complex<double> ans1,ans2,ans_l1,ans_l2;
                        std::tie(ans1,ans2) = find_pair_compl_num(mat, i);
                        std::tie(ans_l1,ans_l2) = find_pair_compl_num(last, i);
                        answer[i] = ans1;
                        answer[i+1] = ans2;
                        if(std::abs(ans1 - ans_l1) > eps || std::abs(ans2 - ans_l2) > eps){
                            check = false;
                            break;
                        }
                        i = i+1;
                        continue;
                    }else{
                        check = false;
                        break;
                    }
                }
            }
            if(check){
                break;
            }
        }
    }
    //std::cout << mat << std::endl;
    return answer;
}


int main(){
    Matrix<double> mat = {{2, -4, 5}, {-5, -2, -3}, {1, -8, -3}};
    //Matrix from example
    //Matrix<double> mat = {{1,3,1},{1,1,4},{4,3,1}};
    std::cout << "Matrix: " << std::endl << mat;
    std::vector<std::complex<double>> ans = find_eigenvalues(mat,0.0001);
    std::cout << "Eigenvalues:" << std::endl;
    for(int i = 0; i < ans.size(); ++i){
        std::cout << "λ" << subs[i+1] << " = ";
        //Normal output complex numbers
        if(ans[i] == 0.){
            std::cout << 0;
        }else{
            bool real = false;
            if(ans[i].real() != 0){
                std::cout << ans[i].real();
                real = true;
            }
            if(ans[i].imag() != 0){
                if(real && ans[i].imag() > 0)
                    std::cout << "+";
                std::cout << ans[i].imag() << "i";
            }
        }
        std::cout << std::endl;
    }
    return 0;
}