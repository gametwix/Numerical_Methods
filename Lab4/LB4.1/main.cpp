#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <iomanip>

double f(double x,double y,double z){
    return -3*y - 2*z/tan(x);
}

double realF(double x){
    return (-0.9783*std::cos(2*x) + 0.4776*std::sin(2*x))/std::sin(x);
}


std::pair<std::vector<double>,std::vector<double>> Eiler(double y0, double z0, double h, double a, double b){
    int n = (b - a) / h;
    double cur_x = a;
    std::vector<double> y(n+1);
    std::vector<double> z(n+1);
    y[0] = y0;
    z[0] = z0;
    for(int i = 1; i <= n; ++i){
        y[i] = y[i - 1] + h*z[i-1];
        z[i] = z[i - 1] + h*f(cur_x,y[i-1],z[i-1]);
        cur_x += h;
    }
    return {y,z};
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

std::pair<std::vector<double>,std::vector<double>> Adams(double y0, double z0, double h, double a, double b){
    std::vector<double> y,z;
    std::tie(y,z) = Runge_Kutta(1,1,h,a,b);
    int n = (b - a) / h;
    std::vector<double> x(4);
    x[0] = a;
    for(int i = 1; i < 4; ++i){
        x[i] = x[i-1]+h;
    }
    for(int i = 4; i <= n; ++i){
        y[i] = y[i - 1] + h/24*(55*z[i-1] - 59*z[i-2] + 37*z[i-3] - 9*z[i-4]);
        z[i] = z[i - 1] + h/24*(55*f(x[3],y[i-1],z[i-1]) - 59*f(x[2],y[i-2],z[i-2]) + 37*f(x[1],y[i-3],z[i-3]) - 9*f(x[0],y[i-4],z[i-4]));
        for(int i = 0; i < 4; ++i){
            x[i] = x[i]+h;
        }
    }
    return {y,z};
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
    double a = 1.;
    double b = 2.;
    std::vector<double> y1,z1;
    std::vector<double> y2,z2;

    std::tie(y1,z1) = Eiler(1,1,h,a,b);    
    std::tie(y2,z2) = Eiler(1,1,h/2,a,b);

    std::cout << std::setprecision(4);
    std::cout << "Eiler:" << std::endl;
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


    std::tie(y1,z1) = Runge_Kutta(1,1,h,a,b);    
    std::tie(y2,z2) = Runge_Kutta(1,1,h/2,a,b);
    std::cout << "Runge_Kutta:" << std::endl;
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


    std::tie(y1,z1) = Adams(1,1,h,a,b);    
    std::tie(y2,z2) = Adams(1,1,h/2,a,b);
    std::cout << "Adams:" << std::endl;
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