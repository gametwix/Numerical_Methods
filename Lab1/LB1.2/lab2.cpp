#include <iostream>
#include <vector>
#include "../matrix.h"

template<class T>
double Tridioiganal_req(Matrix<T> &matrix,Vector<T> &vec,Vector<T> &sol, int step,double prev_p,double prev_q){
    if(step == vec.get_size()){
        return 0;
    }
    double a,b,c,d;
    if(step == 0){
        a = 0;
        b = matrix(0, 0);
        c = matrix(0, 1);
        d = vec(0);
    }else if(step == (vec.get_size() - 1)){
        a = matrix(step, step - 1);
        b = matrix(step, step);
        c = 0;
        d = vec(step);
    }else{
        a = matrix(step, step - 1);
        b = matrix(step, step); 
        c = matrix(step, step + 1);
        d = vec(step);
    }
    double p = ((-1.0)*c)/(b + a*prev_p);
    double q = (d - a*prev_q)/(b+a*prev_p);
    sol(step) = p*Tridioiganal_req(matrix,vec,sol,step+1,p,q) + q;
    return sol(step);
}

template<class T>
Vector<T> Tridioiganal(Matrix<T> &matrix,Vector<T> &b){
    Vector<T> x(b.get_size());
    Tridioiganal_req(matrix,b,x,0,0,0);
    return x;
}

int main(){
    Matrix<double> matrix = {   {-1,-1,0,0,0},
                                {1,-8,1,0,0},
                                {0,-2,-11,5,0},
                                {0,0,3,-14,7},
                                {0,0,0,8,10}};
    
    Vector<double> vec = {-114,81,-8,-38,114};
    Vector<double> sol = Tridioiganal(matrix,vec);
    std::cout << "Matrix: " << std::endl << matrix;
    std::cout << "Vector: " << std::endl << vec;
    std::cout << "Solve: " << std::endl << sol;
    return 0;
}