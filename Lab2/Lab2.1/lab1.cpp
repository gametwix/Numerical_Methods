#include <iostream>
#include <cmath>

double func(double x){
    return x*x*x - 2*x*x - 10*x + 15;
}

double dfunc(double x){
    return 3*x*x - 4*x - 10;
}

double ddfunc(double x){
    return 6*x - 4*x;
}

double phi1(double x){
    return (x*x*x - 2*x*x + 15)/10;
}

double dphi1(double x){
    return (3*x*x - 4*x)/10;
}

double phi2(double x){
    return 2+10/x-15/(x*x);
}

double dphi2(double x){
    return -10/(x*x)+30/(x*x*x);
}


double fixed_point_iteration(double a,double b,double (*function)(double),double (*dfunction)(double),double eps = 0.001){
    double x0 = (b + a) / 2;
    double q = std::max(std::abs(dfunction(a)),std::abs(dfunction(b)));
    double x = function(x0);
    while(q/(1-q)*std::abs(x - x0) > eps){
        x0 = x;
        x = function(x0);
    }
    return x;
}

double newton_method(double a,double b,double (*function)(double),double (*dfunction)(double),double eps = 0.001){
    double x0 = (b + a) / 2;
    double x = x0 - function(x0)/dfunction(x0);
    while(std::abs(x - x0) > eps){
        x0 = x;
        x = x0 - function(x0)/dfunction(x0);
    }
    return x;
}

int main(){
    std::cout << "Функция: " << std::endl;
    std::cout << "f(x) = x^3 - 2x^2 - 10x + 15 = 0" << std::endl;
    std::cout << "Метод простых итераций: " << std::endl;
    std::cout << "x1 = " << fixed_point_iteration(1,2,phi1,dphi1) << std::endl;
    std::cout << "x2 = " << fixed_point_iteration(3,4,phi2,dphi2) << std::endl;
    std::cout << "Метод Ньютона: " << std::endl;
    std::cout << "x1 = " << newton_method(0,2,func,dfunc) << std::endl;
    std::cout << "x2 = " << newton_method(3,6,func,dfunc) << std::endl;
    return 0;
}