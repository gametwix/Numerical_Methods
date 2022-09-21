#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

double func_1(double x1, double x2){
    return x1*x1/9 + x2*x2*4/9 - 1;
}

double func_2(double x1, double x2){
    return 3*x2 - std::exp(x1) - x1;
}

double dfunc_1_1(double x1, double x2){
    return 2*x1/(3*3);
}

double dfunc_1_2(double x1, double x2){
    return x2*8/9;
}

double dfunc_2_1(double x1, double x2){
    return -std::exp(x1) - 1;
}

double dfunc_2_2(double x1, double x2){
    return 3;
}


double phi1(double x1, double x2,double a, double b){
    return x1 - a*func_1(x1,x2) - b*func_2(x1,x2);
}

double dphi1_1(double x1,double x2,double a, double b, double eps = 0.001){
    return (phi1(x1+eps,x2,a,b)-phi1(x1,x2,a,b))/eps;
}

double dphi1_2(double x1,double x2,double a, double b, double eps = 0.001){
    return (phi1(x1,x2+eps,a,b)-phi1(x1,x2,a,b))/eps;
}

double phi2(double x1, double x2,double a, double b){
    return x2 - a*func_1(x1,x2) - b*func_2(x1,x2);
}

double dphi2_1(double x1,double x2,double a, double b, double eps = 0.001){
    return (phi2(x1+eps,x2,a,b)-phi2(x1,x2,a,b))/eps;
}

double dphi2_2(double x1,double x2,double a, double b, double eps = 0.001){
    return (phi2(x1,x2+eps,a,b)-phi2(x1,x2,a,b))/eps;
}

double max_jphi_sum(double x1,double x2,double a, double b,double c, double d){
    return std::max({std::abs(dphi1_1(x1,x2,a,b))+std::abs(dphi1_2(x1,x2,a,b)),
                    std::abs(dphi2_1(x1,x2,c,d))+std::abs(dphi2_2(x1,x2,c,d))});
}


std::vector<double> fixed_point_iteration(double a1,double b1,double a2,double b2, double eps = 0.001){
    double x1_0 = (a1+b1)/2, x2_0= (a2+b2)/2;
    double j1_1 = dfunc_1_1(x1_0,x2_0);
    double j1_2 = dfunc_1_2(x1_0,x2_0);
    double j2_1 = dfunc_2_1(x1_0,x2_0);
    double j2_2 = dfunc_2_2(x1_0,x2_0);
    double det = j1_1*j2_2 - j2_1*j1_2;
    double tmp = j1_1;
    j1_1 = j2_2 / det;
    j2_2 = tmp / det;
    j2_1 /= -det;
    j1_2 /= -det;
    double x1, x2;
    x1 = phi1(x1_0,x2_0,j1_1,j1_2);
    x2 = phi2(x1_0,x2_0,j2_1,j2_2);
    double q = std::max({max_jphi_sum(a1,a2,j1_1,j1_2,j2_1,j2_2),
                        max_jphi_sum(a1,b2,j1_1,j1_2,j2_1,j2_2),
                        max_jphi_sum(b1,a2,j1_1,j1_2,j2_1,j2_2),
                        max_jphi_sum(b1,b2,j1_1,j1_2,j2_1,j2_2)});
    std::cout <<"q=" << q << std::endl;
    while(q/(1-q)*std::max(std::abs(x1 - x1_0),std::abs(x2 - x2_0)) > eps){
        x1_0 = x1;
        x2_0 = x2;
        x1 = phi1(x1_0,x2_0,j1_1,j1_2);
        x2 = phi2(x1_0,x2_0,j2_1,j2_2);
    }
    return {x1,x2};
}


double Det2_2(std::vector<std::vector<double>> vec){
    if(vec.size() == 2 && vec[0].size() == 2){
        return vec[0][0]*vec[1][1] - vec[0][1]*vec[1][0];
    }else{
        throw std::out_of_range("arguments out of range");
    }
}





std::vector<double> newton_method(double a1,double b1,double a2,double b2,double eps = 0.001){
    double x1_0 = (a1+b1)/2, x2_0= (a2+b2)/2;
    double x1, x2;
    std::vector<std::vector<double>> J,A1,A2;
    J = {{dfunc_1_1(x1_0,x2_0),dfunc_1_2(x1_0,x2_0)},{dfunc_2_1(x1_0,x2_0),dfunc_2_2(x1_0,x2_0)}};
    A1 = {{func_1(x1_0,x2_0),dfunc_1_2(x1_0,x2_0)},{func_2(x1_0,x2_0),dfunc_2_2(x1_0,x2_0)}};
    A2 = {{dfunc_1_1(x1_0,x2_0),func_1(x1_0,x2_0)},{dfunc_2_1(x1_0,x2_0),func_2(x1_0,x2_0)}};
    x1 = x1_0 - Det2_2(A1)/Det2_2(J);
    x2 = x2_0 - Det2_2(A2)/Det2_2(J);

    while(std::max(std::abs(x1 - x1_0),std::abs(x2 - x2_0)) > eps){
        x1_0 = x1;
        x2_0 = x2;
        J = {{dfunc_1_1(x1_0,x2_0),dfunc_1_2(x1_0,x2_0)},{dfunc_2_1(x1_0,x2_0),dfunc_2_2(x1_0,x2_0)}};
        A1 = {{func_1(x1_0,x2_0),dfunc_1_2(x1_0,x2_0)},{func_2(x1_0,x2_0),dfunc_2_2(x1_0,x2_0)}};
        A2 = {{dfunc_1_1(x1_0,x2_0),func_1(x1_0,x2_0)},{dfunc_2_1(x1_0,x2_0),func_2(x1_0,x2_0)}};
        x1 = x1_0 - Det2_2(A1)/Det2_2(J);
        x2 = x2_0 - Det2_2(A2)/Det2_2(J);
    }
    return {x1,x2};
}

int main(){
    std::cout << "Функция: " << std::endl;
    std::cout << "f_1(x) = x_1^2/9 + 4*x_2^2/9 - 1 = 0" << std::endl;
    std::cout << "f_2(x) = 3*x_2 - e^x_1 - x_1 = 0" << std::endl;
    std::cout << "Метод простых итераций: " << std::endl;
    std::vector<double> itr = fixed_point_iteration(1,2,1,2);
    for(int i = 0; i < 2; ++i){
        std::cout << "x_" << i+1 << " = " << itr[i] << std::endl;
    }
    std::cout << "Метод Ньютона: " << std::endl;
    std::vector<double> newton = newton_method(0,1.5,0,1.5);
    for(int i = 0; i < 2; ++i){
        std::cout << "x_" << i+1 << " = " << newton[i] << std::endl;
    }
    return 0;
}