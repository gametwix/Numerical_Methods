#include <iostream>
#include <vector>
#include "square_matrix.h"

int main(){
    std::vector<std::vector<double>> A{ {-1,-3,-4,0},
                                        {3,7,-8,3},
                                        {1,-6,2,5},
                                        {-8,-4,-1,-1}};
    SquareMatrix matrix(A);
    std::vector<double> b{-3,30,-90,12};
    std::cout << "Matrix and vector: " << std::endl;
    for(int i = 0; i < matrix.get_size(); ++i){
        for(int j = 0; j < matrix.get_size(); ++j){
            std::cout << matrix.get_elem(i,j) << "\t";
        }
        std::cout <<"= "<<b[i]<< std::endl;
    }
    std::cout << "L: " << std::endl;
    for(int i = 0; i < matrix.get_size(); ++i){
        for(int j = 0; j < matrix.get_size(); ++j){
            std::cout << matrix.get_L_elem(i,j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "U: " << std::endl;
    for(int i = 0; i < matrix.get_size(); ++i){
        for(int j = 0; j < matrix.get_size(); ++j){
            std::cout << matrix.get_U_elem(i,j) << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "Det: " << matrix.get_det() << std::endl;

    
    
    std::vector<double> x = matrix.find_root(b);
    std::cout << "Root: "  << std::endl;
    for(int i = 0; i < x.size();++i){
        std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout <<"Inverse:"<< std::endl;
    SquareMatrix bmat(matrix.get_inverse());
    for(int i = 0; i < bmat.get_size(); ++i){
        for(int j = 0; j < bmat.get_size(); ++j){
            std::cout << bmat.get_elem(i,j) << " \t";
        }
        std::cout << std::endl;
    }
    
    return 0;
}