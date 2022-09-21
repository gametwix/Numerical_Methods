#include <iostream>
#include <vector>
#include <cmath>

double F(double x){
    return std::tan(x) + x;
    //return std::log(x);
    //return std::sin(x*M_PI/6);
}

double Lagrange(double X,const std::vector<double> &Xi,const std::vector<double> &Yi){
    double sum = 0;
    for(int i = 0; i < Xi.size();++i){
        double mult = 1;
        for(int j = 0; j < Xi.size();++j){
            if(i == j) continue;
            mult *= (X - Xi[j])/(Xi[i] - Xi[j]);
        }
        sum += Yi[i]*mult;
    }
    return sum;
}

double DivDiff(const std::vector<double> &Xi,const std::vector<double> &Yi,int i,int j){
    if(i == j){
        return Yi[i];
    }else{
        return (DivDiff(Xi,Yi,i,j-1) - DivDiff(Xi,Yi,i+1,j))/(Xi[i] - Xi[j]);
    }
}


double Newton(double X,const std::vector<double> &Xi,const std::vector<double> &Yi){
    double sum = 0;
    double mult = 1;
    for(int i = 0; i < Xi.size();++i){
        sum += mult*DivDiff(Xi,Yi,0,i);
        mult *= X - Xi[i];
    }
    return sum;
}
int main(){
    std::vector<double> Xi1 = {0,M_PI/8,2*M_PI/8,3*M_PI/8};
    std::vector<double> Yi1 = {F(0),F(M_PI/8),F(2*M_PI/8),F(3*M_PI/8)};
    std::vector<double> Xi2 = {0,M_PI/8,M_PI/3,3*M_PI/8};
    std::vector<double> Yi2 = {F(0),F(M_PI/8),F(M_PI/3),F(3*M_PI/8)};
    double X = 3*M_PI/16;

    std::cout << "1) ";
    for(int i = 0; i< Xi1.size();++i){
        std::cout << (i == 0 ? "" : ", ") << Xi1[i];
    }
    std::cout << std::endl;

    std::cout << "2) ";
    for(int i = 0; i< Xi2.size();++i){
        std::cout << (i == 0 ? "" : ", ") << Xi2[i];
    }
    std::cout << std::endl;

    std::cout << "X* = "<< X <<std::endl;

    std::cout << "1)" << std::endl;
    std::cout << "Lagrange: Y* = " <<  Lagrange(X,Xi1,Yi1) << " Eps = " << std::abs(Lagrange(X,Xi1,Yi1) - F(X)) << std::endl;
    std::cout << "Newton: Y* = " << Newton(X,Xi1,Yi1) << " Eps = " << std::abs(Newton(X,Xi1,Yi1) - F(X))  << std::endl;
    std::cout << "2)" << std::endl;
    std::cout << "Lagrange: Y* = "<< Lagrange(X,Xi2,Yi2) << " Eps = " << std::abs(Lagrange(X,Xi2,Yi2) - F(X))  << std::endl;
    std::cout << "Newton: Y* = " << Newton(X,Xi2,Yi2) << " Eps = " << std::abs(Newton(X,Xi2,Yi2) - F(X))  << std::endl;
    return 0;
}