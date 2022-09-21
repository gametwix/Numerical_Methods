#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>
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

double f(double x, double y,double z){
    return (2*z + std::exp(x)*y) / (std::exp(x) + 1);
}

double p(double x){
    return -2 / (std::exp(x) + 1);
}

double q(double x){
    return (-std::exp(x)) / (std::exp(x) + 1);
}

double realF(double x){
    return std::exp(x) - 1;
}

std::pair<double,double> delt(double x, double y, double z, double h){
    std::vector<double> K(4);
    std::vector<double> L(4);
    K[0] = h*z;
    L[0] = h*f(x,y,z);
    for(int i = 1; i < 4; ++i){
        K[i] = h*(z + L[i-1]/2);
        L[i] = h*f(x + h/2, y + K[i-1]/2,z + L[i-1]/2);
    }
    return {(K[0] + 2*K[1] + 2*K[2] + K[3])/6, (L[0] + 2*L[1] + 2*L[2] + L[3])/6};
}



std::pair<std::vector<double>,std::vector<double>> Runge_Kutta(double y0, double z0, double h, double a, double b){
    int n = (b - a) / h;
    double cur_x = a;
    std::vector<double> y(n+1);
    std::vector<double> z(n+1);
    y[0] = y0;
    z[0] = z0;
    for(int i = 1; i <= n; ++i){
        double dy, dz;
        std::tie(dy, dz) = delt(cur_x, y[i-1], z[i-1], h);
        y[i] = y[i - 1] + dy;
        z[i] = z[i - 1] + dz;
        cur_x += h;
    }
    return {y,z};
}

double F(double ay, double by, double n, double a, double b){
    std::vector<double> y,z;
    std::tie(y,z) = Runge_Kutta(ay,n,0.1,a,b);
    return y[y.size() - 1] - by;
}

std::pair<std::vector<double>,std::vector<double>> shooting(double a, double b, double ay, double by,double h, double eps){
    std::vector<double> nt(3);
    std::vector<double> nF(3);
    std::vector<double> y,z;
    nt[1] = 0.1;
    nt[2] = 3;
    nF[1] = F(ay,by,nt[1],a,b);
    nF[2] = F(ay,by,nt[2],a,b);
    do{
        nt[0] = nt[1];
        nt[1] = nt[2];
        nF[0] = nF[1];
        nF[1] = nF[2];
        nt[2] = nt[1] -  (nt[1] - nt[0]) / (nF[1] - nF[0])*nF[1];
        nF[2] = F(ay,by,nt[2],a,b);
        //std::cout << nt[2] << " " << nF[2] << " " << std::abs(nF[2])<< std::endl;
    }while(std::abs(nF[2]) >= eps);
    return Runge_Kutta(ay,nt[2],h,a,b);
}


std::vector<double> Raznost(double a, double b, double h, double ay, double by){
    
    int n = std::abs(b - a) / h;
    Matrix<double> to_solv(n-1,n-1);
    Vector<double> vec(n-1);
    std::vector<double> y(n+1);
    y[0] = ay;
    y[n] = by;
    double cut_x = a+h;
    to_solv(0,0) = -2 + h*h*q(cut_x);
    to_solv(0,1) = 1 + h*p(cut_x)/2;
    vec(0) = p(cut_x)*h*ay/2 - ay;

    for(int i = 1; i < n - 2; ++i){
        //std::cout << to_solv << std::endl;
        //std::cout << i << std::endl;
        cut_x = a+h;
        to_solv(i,i-1) = (1 - h*p(cut_x)/2);
        to_solv(i,i) = -2 + h*h*q(cut_x);
        to_solv(i,i+1) = (1 + h*p(cut_x)/2);
        vec(i) = 0;
    }
    //std::cout << to_solv << std::endl;
    cut_x = a+h;
    to_solv(n - 2,n - 3) = (1 - h*p(cut_x)/2);
    //std::cout << to_solv << std::endl;
    to_solv(n - 2,n - 2) = -2 + h*h*q(cut_x);
   
    vec(n - 2) = -by - by*p(cut_x)*h/2;
    Vector<double> ans_c = Tridioiganal(to_solv, vec);
    //std::cout << ans_c << std::endl;

    for(int i = 1; i < y.size() - 1; ++i){
        y[i] = ans_c(i-1);
    }
    return y;
}

double eps(std::vector<double> y, double a, double b, double h){
    double sum = 0;
    int n = (b - a) / h;
    double cur_x = a;
    for(int i = 0; i <= n; ++i){
        sum += std::pow(y[i] - realF(cur_x),2);
        cur_x += h;
    }
    return std::sqrt(sum);
}

double Runge_Romberg(double h1, double h2,std::vector<double> I1,std::vector<double> I2,double p){
    if(h1 > h2){
        double sum = 0;
        for(int i = 0; i < I1.size(); ++i){
            sum += std::pow(I1[i] - I2[2*i],2);
        }
        sum = std::sqrt(sum);
        return sum / (std::pow((h2 / h1),p) - 1);
    }else{
        double sum = 0;
        for(int i = 0; i < I2.size(); ++i){
            sum += std::pow(I1[2*i] - I2[i],2);
        }
        sum = std::sqrt(sum);
        return sum / (std::pow((h1 / h2),p) - 1);
    }
}

int main(){
    double h = 0.1;
    double a = 0, b = 1.;
    double ay = 0, by = 1.72;
    std::vector<double> y1,z1,y2,z2;
   
    std::tie(y1,z1) = shooting(a,b,ay,by,h,0.0001);
    std::tie(y2,z2) = shooting(a,b,ay,by,h/2,0.0001);
    std::cout << std::setprecision(4);
    std::cout << "Метод стрельбы :" << std::endl;
    std::cout << "x:" << std::defaultfloat;
    for(int i = 0; i <= (b - a) / h; ++i){
        std::cout << std::setw(7) << a + i*h;
    }
    std::cout << std::endl;
    std::cout << "y:";
    std::cout << std::fixed;
    for(int i = 0; i < y1.size(); ++i){
        std::cout << std::setw(7) << y1[i];
    }
    std::cout << std::endl;
    std::cout << "Погрешность Рунге-Ромберга: " << std::setw(7) << Runge_Romberg(h,h/2,y1,y2,2)<< std::endl;
    std::cout << "Погрешность относительно точного решения: " << std::setw(7) << eps(y1,a,b,h) << std::endl << std::endl;


    y1 = Raznost(a,b,h,ay,by);    
    y2 = Raznost(a,b,h,ay,by);
    std::cout << "Конечно-разностный метод:" << std::endl;
    std::cout << "x:" << std::defaultfloat;
    for(int i = 0; i <= (b - a) / h; ++i){
        std::cout << std::setw(7) << a + i*h;
    }
    std::cout << std::endl;
    std::cout << "y:";
    std::cout << std::fixed;
    for(int i = 0; i < y1.size(); ++i){
        std::cout << std::setw(7) << y1[i];
    }
    std::cout << std::endl;
    std::cout << "Погрешность Рунге-Ромберга: " << std::setw(7) << Runge_Romberg(h,h/2,y1,y2,2)<< std::endl;
    std::cout << "Погрешность относительно точного решения: " << std::setw(7) <<eps(y1,a,b,h) << std::endl << std::endl;

    return 0;
}