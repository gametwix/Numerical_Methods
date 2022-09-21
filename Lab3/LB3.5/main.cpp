#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

double f(double x){
    return 1 / (std::pow(x,4) + 16);
}

double Runge_Romberg(double h1, double h2,double I1,double I2,double p){
    if(h1 > h2){
        return I1 + (I1 - I2) / (std::pow((h2 / h1),p) - 1);
    }else{
        return I2 + (I2 - I2) / (std::pow((h1 / h2),p) - 1);
    }
}

double Rectangle_method(double X1,double Xk, double h){
    int k = (Xk - X1) / h;
    double sum = 0;
    double LX = X1, NX = X1 + h;
    for(int i = 0; i < k; ++i){
        sum += h*f((LX+NX) / 2);
        LX = NX;
        NX += h;
    }
    return sum;
}

double Trapezoidal_method(double X1,double Xk, double h){
    int k = (Xk - X1) / h;
    double sum = 0;
    double LX = X1, NX = X1 + h;
    for(int i = 0; i < k; ++i){
        sum += h*(f(LX)+f(NX))/2;
        LX = NX;
        NX += h;
    }
    return sum;
}

double Simpson_method(double X1,double Xk, double h){
    int k = (Xk - X1) / (2*h);
    double sum = 0;
    double LX = X1, MX = X1 + h, NX = X1 + 2*h;
    for(int i = 0; i < k; ++i){
        sum += h*(f(LX)+ 4*f(MX) + f(NX))/3;
        LX = NX;
        MX += 2*h;
        NX += 2*h;
    }
    return sum;
}

int main(){
    double X1 = 0.;
    double Xk = 2.;
    double h1 = 0.5;
    double h2 = 0.25;
    std::cout << "h1 = " << h1 << std::endl;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Rectangle_method: "  << Rectangle_method(X1,Xk,h1) << std::endl;
    std::cout << "Trapezoidal_method: "  << Trapezoidal_method(X1,Xk,h1) << std::endl;
    std::cout << "Simpson_method: " << Simpson_method(X1,Xk,h1) << std::endl<< std::endl;

    std::cout << std::defaultfloat << "h2 = " << h2 << std::endl << std::fixed;
    std::cout << "Rectangle_method: " << Rectangle_method(X1,Xk,h2) << std::endl;
    std::cout << "Trapezoidal_method: " << Trapezoidal_method(X1,Xk,h2) << std::endl;
    std::cout << "Simpson_method: "  << Simpson_method(X1,Xk,h2) << std::endl<< std::endl;

    std::cout << "Runge_Romberg" << std::endl<< std::endl;
    std::cout << "Rectangle_methods" << std::endl;
    std::cout <<"Answer: " << Runge_Romberg(h1,h2,Rectangle_method(X1,Xk,h1),Rectangle_method(X1,Xk,h2),2) << std::endl;
    std::cout <<"Eps1: " << std::abs(Runge_Romberg(h1,h2,Rectangle_method(X1,Xk,h1),Rectangle_method(X1,Xk,h2),2) - Rectangle_method(X1,Xk,h1)) <<"\tEps2: "
     << std::abs(Runge_Romberg(h1,h2,Rectangle_method(X1,Xk,h1),Rectangle_method(X1,Xk,h2),2) - Rectangle_method(X1,Xk,h2))<< std::endl<< std::endl;
    std::cout << "Trapezoidal_method" << std::endl;
    std::cout <<"Answer: " <<Runge_Romberg(h1,h2,Trapezoidal_method(X1,Xk,h1),Trapezoidal_method(X1,Xk,h2),2) << std::endl;
    std::cout <<"Eps1: " << std::abs(Runge_Romberg(h1,h2,Trapezoidal_method(X1,Xk,h1),Trapezoidal_method(X1,Xk,h2),2) - Trapezoidal_method(X1,Xk,h1)) <<"\tEps2: "
     << std::abs(Runge_Romberg(h1,h2,Trapezoidal_method(X1,Xk,h1),Trapezoidal_method(X1,Xk,h2),2) - Trapezoidal_method(X1,Xk,h2))<< std::endl<< std::endl;
    std::cout << "Simpson_method" << std::endl;
    std::cout <<"Answer: " <<Runge_Romberg(h1,h2,Simpson_method(X1,Xk,h1),Simpson_method(X1,Xk,h2),2) << std::endl;
    std::cout <<"Eps1: " << std::abs(Runge_Romberg(h1,h2,Simpson_method(X1,Xk,h1),Simpson_method(X1,Xk,h2),2) - Simpson_method(X1,Xk,h1)) <<"\tEps2: "
     << std::abs(Runge_Romberg(h1,h2,Simpson_method(X1,Xk,h1),Simpson_method(X1,Xk,h2),2) - Simpson_method(X1,Xk,h2))<< std::endl;
    return 0;
}