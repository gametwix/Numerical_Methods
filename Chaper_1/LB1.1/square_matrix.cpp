#include "square_matrix.h"

bool check_is_matrix(const std::vector<std::vector<double>> &A){
    if(A.empty()){
        return false;
    }
    int rows = A.size();
    int colums = A[0].size();
    for(int i = 1; i < rows; ++i){
        if(A[i].size() != colums){
            return false;
        }
    }
    return true;
}

bool check_is_square_matrix(const std::vector<std::vector<double>> &A){
    if(!check_is_matrix(A)){
        return false;
    }
    int rows = A.size();
    int colums = A[0].size();
    return rows == colums;
}

SquareMatrix::SquareMatrix():size(0){}

SquareMatrix::SquareMatrix(const int _size):size(_size){
    elements = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
    L = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
    U = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
}

SquareMatrix::SquareMatrix(const std::vector<std::vector<double>> _elements){
    if(!check_is_square_matrix(_elements)){
        std::cerr << "Error: Vector is not a square matrix" << std::endl;
        throw -1;
    }
    size = _elements.size();
    elements = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
    L = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
    U = std::vector<std::vector<double>>(size,std::vector<double>(size,0.0));
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            elements[i][j] = _elements[i][j];
        }
    }
    find_LU();
    
}

SquareMatrix::SquareMatrix(const SquareMatrix &_matrix){
    size = _matrix.size;
    elements = _matrix.elements;
    L = _matrix.L;
    U = _matrix.U;
}

int SquareMatrix::get_size(){return size;}

double SquareMatrix::get_elem(int i, int j){return elements[i][j];}

void SquareMatrix::set_elem(int i,int j,double elem){elements[i][j] = elem;}

double SquareMatrix::get_L_elem(int i, int j){return L[i][j];}
double SquareMatrix::get_U_elem(int i, int j){return U[i][j];}

void SquareMatrix::find_LU(){
    U = elements;
    for(int k = 0; k < size;++k){
        L[k][k] = 1;
        for(int i = k+1;i<size;++i){
            L[i][k] = U[i][k] / U[k][k];
            for(int j = k; j < size; ++j){
                U[i][j] -= U[k][j]*L[i][k];
            }
        }
    }
}

double SquareMatrix::get_det(){
    if(size < 1){
        std::cerr << "Error: Can't find determinant in empty matrix" << std::endl;
        throw -1;
    }
    double ans = 1;
    for(int i = 0; i < size; ++i){
        ans *= U[i][i];
    }
    return ans;
}

std::vector<double> SquareMatrix::find_root(const std::vector<double> B){
    if(size == 0){
        std::cerr << "Error: Can't find root with empty matrix" << std::endl;
        throw -1;
    }
    if(size != B.size()){
        std::cerr << "Error: Different size of matrix and vector" << std::endl;
        throw -1;
    }
    std::vector<double> y(size);

    for(int i = 0;i < size; ++i){
        y[i] = B[i];
        for(int j = 0; j< i;++j){
            y[i] -= y[j]*L[i][j];
        }
        y[i] /= L[i][i];
    }
    std::vector<double> x(size);
    for(int i = size - 1;i >= 0; --i){
        x[i] = y[i];
        for(int j = size - 1; j > i;--j){
            x[i] -= x[j]*U[i][j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

SquareMatrix SquareMatrix::get_inverse(){
    SquareMatrix answer(size);
    for(int i = 0; i < size;++i){
        std::vector<double> b(size,0.0);
        b[i] = 1.0;
        std::vector<double> x = find_root(b);
        for(int j = 0; j < size; ++j){
            answer.set_elem(j,i,x[j]);
        }
    }
    answer.find_LU();
    return answer;
}