#include <iostream>
#include <cmath>
#include "../matrix.h"

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
Matrix<T> find_Q(Matrix<T> mat){
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
    return Q;
}


int main(){
    //Matrix<double> mat = {{2, -4, 5}, {-5, -2, -3}, {1, -8, -3}};
    Matrix<double> mat = {{1,3,1},{1,1,4},{4,3,1}};
    std::cout << find_Q(mat) << std::endl;
    return 0;
}