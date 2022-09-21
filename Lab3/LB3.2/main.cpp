#include <iostream>
#include <cmath>
#include <vector>
#include "../../Lab1/matrix.h"

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

double h(const std::vector<double> Xi, int i){
    return Xi[i] - Xi[i-1];
}

std::vector<double> gen_c(const std::vector<double> &Xi,const std::vector<double> &Fi){
    Matrix<double> to_solve_c(Xi.size() - 2,Xi.size() - 2);
    Vector<double> vec(Xi.size() - 2);
    std::vector<double> c(Xi.size() - 1);

    to_solve_c(0, 2 - 2) = 2*(h(Xi,1) + h(Xi,2));
    to_solve_c(0, 3 - 2) = h(Xi,2);
    vec(0) = 3*((Fi[2] - Fi[1])/ h(Xi,2) - (Fi[1] - Fi[0])/ h(Xi,1));

    for(int i = 3; i < Xi.size() - 1; ++i){
        to_solve_c(i - 2, i - 1 - 2) =  h(Xi, i - 1);
        to_solve_c(i - 2, i - 2) = 2*(h(Xi, i - 1) + h(Xi, i));
        to_solve_c(i - 2, i + 1 - 2) = h(Xi, i);
        vec(i - 2) = 3*((Fi[i] - Fi[i - 1]) / h(Xi, i) - (Fi[i - 1] - Fi[i - 2]) / h(Xi, i - 1));
    }

    to_solve_c(Xi.size() - 1 - 2, Xi.size() - 1 - 1 - 2) =  h(Xi, Xi.size() - 1 - 1);
    to_solve_c(Xi.size() - 1 - 2, Xi.size() - 1 - 2) = 2*(h(Xi, Xi.size() - 1 - 1) + h(Xi, Xi.size() - 1));
    vec(Xi.size() - 1 - 2) = 3*((Fi[Xi.size() - 1] - Fi[Xi.size() - 1 - 1]) / h(Xi, Xi.size() - 1) - (Fi[Xi.size() - 1 - 1] - Fi[Xi.size() - 1 - 2]) / h(Xi, Xi.size() - 1 - 1));
    Vector<double> ans_c = Tridioiganal(to_solve_c, vec);

    for(int i = 1; i < Xi.size() - 1; ++i){
        c[i] = ans_c(i - 1);
    }
    return c;
}

double S(double X, const std::vector<double> &Xi,
                   const std::vector<double> &a,
                   const std::vector<double> &b,
                   const std::vector<double> &c,
                   const std::vector<double> &d){
    int m = 0;
    double Xm =0;
    for(int i = 0; i < Xi.size() - 1; ++i){
        if(Xi[i] <= X && X <= Xi[i+1]){
            m = i;
            Xm = Xi[i];
            break;
        }
    }
    return a[m] + b[m]*(X - Xm) + c[m]*(X - Xm)*(X - Xm) + d[m]*(X - Xm)*(X - Xm)*(X - Xm);
}


int main(){
    std::vector<double> Xi = {0.0, 0.9, 1.8, 2.7, 3.6};
    std::vector<double> Fi = {0.0, 0.72235, 1.5609, 2.8459, 7.7275};
    //std::vector<double> Xi = {0.0, 1.0, 2.0, 3.0, 4.0};
    //std::vector<double> Fi = {0.0, 1.8415, 2.9093, 3.1411, 3.2432};
    double X = 1.5;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << " i ";
    for(size_t i = 0; i < 5; ++i){
        std::cout << "   " << i << "   ";
    }
    std::cout << std::endl;
    std::cout << " xᵢ";
    for(size_t i = 0; i < 5; ++i){
        std::cout << std::setw(7) << Xi[i];
    }
    std::cout << std::endl;
    std::cout << " fᵢ";
    for(size_t i = 0; i < 5; ++i){
        std::cout << std::setw(7) << Fi[i];
    }
    std::cout << std::endl << std::endl << "X* = " << std::defaultfloat << X <<  std::endl << std::endl << std::fixed;

    std::vector<double> a(Xi.size() - 1);
    std::vector<double> b(Xi.size() - 1);
    std::vector<double> c(Xi.size() - 1);
    std::vector<double> d(Xi.size() - 1);
    c = gen_c(Xi,Fi);

    for(int i = 0; i < c.size() - 1; ++i){
        a[i] = Fi[i];
        b[i] = (Fi[i+1] - Fi[i]) / h(Xi,i+1) - (c[i+1] + 2*c[i])*h(Xi,i+1)/3;
        d[i] = (c[i+1] - c[i]) / (3*h(Xi,i+1));
    }

    a[c.size() - 1] = Fi[c.size() - 1];
    b[c.size() - 1] = (Fi[c.size()] - Fi[c.size() - 1]) / h(Xi,c.size()) - 2*h(Xi,c.size())*c[c.size() - 1]/3;
    d[c.size() - 1] = (- c[c.size() - 1]) / (3*h(Xi,c.size()));
    

    std::cout << " i    aᵢ     bᵢ     cᵢ     dᵢ  " << std::endl;
    for(int i = 0; i < 4;++i){
        std::cout << " " << i+1 << " "<< std::setw(7) << a[i] << std::setw(7) << b[i] << std::setw(7) << c[i] << std::setw(7) << d[i] << std::endl;
    }

    std::cout << std::endl << "f(x) = " << S(X,Xi,a,b,c,d) << std::endl;
    return 0;
}