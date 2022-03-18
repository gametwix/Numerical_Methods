#include <iostream>
#include <vector>
#include <cmath>

std::vector<double>multiplication(std::vector<std::vector<double>> &matrix,std::vector<double> &b){
    std::vector<double> res(matrix.size());
    for(int i = 0; i < matrix.size();++i){
        res[i] = 0;
        for(int j = 0; j < matrix.size();++j){
            res[i] += matrix[i][j]*b[j];
        }
    }
    return res;
}

std::vector<double> iter(std::vector<std::vector<double>> &matrix,std::vector<double> &b,double eps){
    std::vector<std::vector<double>> alf(matrix.size(),std::vector<double>(matrix.size()));
    std::vector<double> bet(matrix.size());
    for(int i = 0; i < matrix.size();++i){
        for(int j = 0; j < matrix.size();++j){
            if(j == i){
                alf[i][j] = 0;
            }else {
                alf[i][j] = (-1) * matrix[i][j] / matrix[i][i];
            }
        }
        bet[i] = b[i] / matrix[i][i];
    }

    std::vector<double> last_x = bet;
    std::vector<double> cur_x(matrix.size());
    double cur_eps;
    int i = 0;
    do{
        cur_x = multiplication(alf,last_x);
        for(int i = 0; i < matrix.size();++i){
            cur_x[i] += bet[i];
        }
        cur_eps = 0;
        for(int i = 0; i < matrix.size();++i){
            cur_eps += (cur_x[i] - last_x[i])*(cur_x[i] - last_x[i]);
        }
        cur_eps = std::sqrt(cur_eps);
        last_x = cur_x;
        
        ++i;
    }while(cur_eps > eps);
    std::cout << "Count iter:" << i << std::endl;
    return cur_x;
}

int main(){
    std::vector<std::vector<double>> matrix = {{-22,-2,-6,6},{3,-17,-3,7},{2,6,-17,5},{-1,-8,8,23}};
    std::vector<double> b = {96,-26,35,-234};
    double eps = 0.00001;
    std::cout << "Matrix: " << std::endl;
    for(auto &vec : matrix){
        for(auto &elem : vec){
            std::cout << elem << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "B: " << std::endl;
    for(auto &elem : b){
        std::cout << elem << "\t";
    }
    std::cout << std::endl;
    std::cout << "Eps: " << eps << std::endl;
    std::vector<double> x = iter(matrix,b,eps);
    for(int i = 0;i < x.size();++i){
        std::cout << x[i] << "\t";
    }
    std::cout << std::endl;
    return 0;
}