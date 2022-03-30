#include <iostream>
#include <vector>
#include <cmath>
#include "../matrix.h"

template<class T>
Vector<T> zeidel_multiplication(Matrix<T> &matrix,Vector<T> &b, Vector<T> &to_sum){
    Vector<T> res(b.get_size());
    for(int i = 0; i < b.get_size();++i){
        res(i) = b(i);
    }
    for(int i = 0; i < matrix.get_row();++i){
        res(i) = to_sum(i);
        for(int j = 0; j < matrix.get_col();++j){
            res(i) += matrix(i, j)*res(j);
        }
    }
    return res;
}

template<class T>
Vector<T> iter(Matrix<T> &matrix,Vector<T> &b,double eps){
    Matrix<T> alf(matrix.get_row(),matrix.get_col());
    Vector<T> bet(matrix.get_row());
    for(int i = 0; i < matrix.get_row();++i){
        for(int j = 0; j < matrix.get_col();++j){
            if(j == i){
                alf(i, j) = 0;
            }else {
                alf(i, j) = (-1) * matrix(i, j) / matrix(i, i);
            }
        }
        bet(i) = b(i) / matrix(i, i);
    }

    Vector<T> last_x = bet;
    Vector<T> cur_x(matrix.get_row());
    double cur_eps;
    int i = 0;
    do{
        cur_x = alf * last_x;
        for(int i = 0; i < matrix.get_row();++i){
            cur_x(i) += bet(i);
        }
        cur_eps = 0;
        for(int i = 0; i < matrix.get_row();++i){
            cur_eps += (cur_x(i) - last_x(i))*(cur_x(i) - last_x(i));
        }
        cur_eps = std::sqrt(cur_eps);
        last_x = cur_x;
        ++i;
    }while(cur_eps > eps);
    std::cout << "Count iter:" << i << std::endl;
    return cur_x;
}

template<class T>
Vector<T> zeidel(Matrix<T> &matrix,Vector<T> &b,double eps){
    Matrix<T> alf(matrix.get_row(),matrix.get_col());
    Vector<T> bet(matrix.get_row());
    for(int i = 0; i < matrix.get_row();++i){
        for(int j = 0; j < matrix.get_col();++j){
            if(j == i){
                alf(i, j) = 0;
            }else {
                alf(i, j) = (-1) * matrix(i, j) / matrix(i, i);
            }
        }
        bet(i) = b(i) / matrix(i, i);
    }

    Vector<T> last_x = bet;
    Vector<T> cur_x(matrix.get_row());
    double cur_eps;
    int i = 0;
    do{
        cur_x = zeidel_multiplication(alf,last_x,bet);
        cur_eps = 0;
        for(int i = 0; i < matrix.get_row();++i){
            cur_eps += (cur_x(i) - last_x(i))*(cur_x(i) - last_x(i));
        }
        cur_eps = std::sqrt(cur_eps);
        last_x = cur_x;
        ++i;
    }while(cur_eps > eps);
    std::cout << "Count iter:" << i << std::endl;
    return cur_x;
}

int main(){
    Matrix<double> matrix ={ {-22,-2,-6,6},
                             {3,-17,-3,7},
                             {2,6,-17,5},
                             {-1,-8,8,23}};
    Vector<double> vec = {96,-26,35,-234};
    double eps = 0.001;
    std::cout << "Matrix: " << std::endl << matrix;
    std::cout << "Vector: " << std::endl << vec;
    std::cout << "Eps: " << eps << std::endl;
    std::cout << "Iterations: " << std::endl;
    Vector<double> ans = iter(matrix,vec,eps);
    std::cout << "Answer: " << std::endl << ans;
    std::cout << "Zeidel: " << std::endl;
    ans = zeidel(matrix,vec,eps);
    std::cout << "Answer: " << std::endl << ans;
    return 0;
}