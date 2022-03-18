#include <iostream>
#include <vector>
#include "square_matrix.h"

int main(){
    std::vector<std::vector<double>> A{ {-1,-3,-4,0},
                                        {3,7,-8,3},
                                        {1,-6,2,5},
                                        {-8,-4,-1,-1}};
    std::vector<double> b{-3,30,-90,12};
    SquareMatrix matrix(4);
    Vector<double> vec(4);

    for(int i = 0; i < 4; ++i){
        for(int j = 0; j < 4; ++j){
            matrix(i,j) = A[i][j];
        }
        vec(i) = b[i];
    }
    matrix.find_LU();
    std::cout << "Matrix and vector: " << std::endl;
    std::cout << matrix;
    std::cout << "L: " << std::endl;
    std::cout << matrix.L;
    std::cout << "U: " << std::endl;
    std::cout << matrix.U;
    std::cout << "Det: " << matrix.get_det() << std::endl;
    Vector<double> x = matrix.find_root(vec);
    std::cout << "Root: "  << std::endl;
    std::cout << x;
    std::cout << std::endl;
    std::cout <<"Inverse:"<< std::endl;
    SquareMatrix bmat(matrix.get_inverse());
    std::cout << bmat;
    
    return 0;
}