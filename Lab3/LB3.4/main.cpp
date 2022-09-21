#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

double deriv_1_1(double X, const std::vector<double> &Xi, const std::vector<double> &Yi,int i){
    return (Yi[i+1] - Yi[i]) / (Xi[i+1] - Xi[i]);
}

double deriv_1_2(double X, const std::vector<double> &Xi, const std::vector<double> &Yi,int i){
    return (Yi[i+1] - Yi[i]) / (Xi[i+1] - Xi[i]) + ( (Yi[i+2] - Yi[i+1]) / (Xi[i+2] - Xi[i+1]) - (Yi[i+1] - Yi[i]) / (Xi[i+1] - Xi[i]))*(2*X - Xi[i] - Xi[i+1])/(Xi[i+2] - Xi[i]);
}

double deriv_2(double X, const std::vector<double> &Xi, const std::vector<double> &Yi,int i){
    return 2*( (Yi[i+2] - Yi[i+1]) / (Xi[i+2] - Xi[i+1]) - (Yi[i+1] - Yi[i]) / (Xi[i+1] - Xi[i]))/(Xi[i+2] - Xi[i]);
}

int main(){
    std::vector<double> Xi = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> Yi = {1.0, 2.6931, 4.0986, 5.3863, 6.6094};
    double X = 3.0;

    //Print input data
    std::cout << " i ";
    std::cout << std::fixed << std::setprecision(4);
    for(size_t i = 0; i < Xi.size(); ++i){
        std::cout << "   " << i << "   ";
    }
    std::cout << std::endl;
    std::cout << " xᵢ";
    for(size_t i = 0; i < Xi.size(); ++i){
        std::cout << std::setw(7) << Xi[i];
    }
    std::cout << std::endl;
    std::cout << " yᵢ";
    for(size_t i = 0; i < Yi.size(); ++i){
        std::cout << std::setw(7) << Yi[i];
    }
    std::cout << std::endl << "X* = " << std::defaultfloat << X << std::fixed << std::endl << std::endl;

    std::cout << "Первые производные с первым порядком точности: " << std::endl;
    for(int i = 0; i < Xi.size() - 1; ++i){
        if(Xi[i] <= X && X <= Xi[i+1]){
            std::cout <<"[" << std::defaultfloat << Xi[i] << "," << Xi[i+1] << "] : y'(x) = " << std::fixed << deriv_1_1(X,Xi,Yi,i) << std::endl;
        }
    }
    std::cout << "Первые производные со вторым порядком точности: " << std::endl;
    for(int i = 0; i < Xi.size() - 1; ++i){
        if(Xi[i] <= X && X <= Xi[i+1]){
            std::cout <<"[" << std::defaultfloat << Xi[i] << "," << Xi[i+1] << "] : y'(x) = " << std::fixed << deriv_1_2(X,Xi,Yi,i) << std::endl;
        }
    }
    std::cout << "Вторые производные: " << std::endl;
    for(int i = 0; i < Xi.size() - 2; ++i){
        if(Xi[i] <= X && X <= Xi[i+1]){
            std::cout <<"[" << std::defaultfloat << Xi[i] << "," << Xi[i+1] << "] : y'(x) = " << std::fixed << deriv_2(X,Xi,Yi,i) << std::endl;
        }
    }

    return 0;
}