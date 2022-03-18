#include <iostream>
#include <vector>

double Tridioiganal_req(std::vector<std::vector<double>> &matrix,std::vector<double> &r,std::vector<double> &x, int step,double prev_p,double prev_q){
    if(step == r.size()){
        return 0;
    }
    double a,b,c,d;
    if(step == 0){
        a = 0;
        b = matrix[0][0];
        c = matrix[0][1];
        d = r[0];
    }else if(step == (r.size() - 1)){
        a = matrix[step][step - 1];
        b = matrix[step][step];
        c = 0;
        d = r[step];
    }else{
        a = matrix[step][step - 1];
        b = matrix[step][step]; 
        c = matrix[step][step + 1];
        d = r[step];
    }
    double p = ((-1.0)*c)/(b + a*prev_p);
    double q = (d - a*prev_q)/(b+a*prev_p);
    x[step] = p*Tridioiganal_req(matrix,r,x,step+1,p,q) + q;
    return x[step];
}

std::vector<double> Tridioiganal(std::vector<std::vector<double>> &matrix,std::vector<double> &b){
    std::vector<double> x(b.size());
    Tridioiganal_req(matrix,b,x,0,0,0);
    return x;
}

int main(){
    std::vector<std::vector<double>> matrix = {{-1,-1,0,0,0},{1,-8,1,0,0},{0,-2,-11,5,0},{0,0,3,-14,7},{0,0,0,8,10}};
    std::vector<double> b = {-114,81,-8,-38,114};
    std::vector<double> x = Tridioiganal(matrix,b);
    std::cout << "Matrix: " << std::endl;
    for(int i = 0; i < matrix.size(); ++i){
        for(int j = 0; j < matrix[i].size(); ++j){
            std::cout << matrix[i][j] << "\t" ;
        }
        std::cout << std::endl;
    }
    std::cout << "b - vector: " << std::endl;
    for(int i = 0; i < b.size();++i){
        std::cout << b[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Solve: " << std::endl;
    for(int i = 0; i < x.size();++i){
        std::cout << x[i] << "\t";
    }
    std::cout << std::endl;

    return 0;
}