#include <iostream>
#include <vector>
#include "square_matrix.h"

int main(){
    SquareMatrix matrix{ {-1.0,-3.0,-4.0,0.0},
                                        {3.0,7.0,-8.0,3.0},
                                        {1.0,-6.0,2.0,5.0},
                                        {-8.0,-4.0,-1.0,-1.0}};                                      
    Vector<double> vec{-3,30,-90,12};
    matrix.find_LU();
    
    std::cout << "Matrix: " << matrix << "Vector: "<< std::endl << vec;
    std::cout << "L: " << std::endl << matrix.L;
    std::cout << "U: " << std::endl << matrix.U;
    std::cout << "Det: " << matrix.get_det() << std::endl;
    Vector<double> x = matrix.find_root(vec);
    std::cout << "Root: " << std::endl << x;
    SquareMatrix bmat(matrix.get_inverse());
    std::cout <<"Inverse:"<< std::endl << bmat;
    
    return 0;
}