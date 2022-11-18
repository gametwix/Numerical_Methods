#pragma once
#include "matrix.h"
#include <iostream>
#include <iomanip>
#include <vector>

class SquareMatrix: public Matrix<double>{
public:
    SquareMatrix(const size_t inp_size);
    SquareMatrix();
    SquareMatrix(std::initializer_list<std::initializer_list<double>> elems):   Matrix(elems),
                                                                                size(elems.size()),
                                                                                U(elems.size(),elems.size()), 
                                                                                L(elems.size(),elems.size()){}

    size_t get_size() const;

    double get_L_elem(int i, int j) const;
    double get_U_elem(int i, int j) const;

    double get_det();
    SquareMatrix get_inverse();
    Vector<double> find_root(Vector<double> B);  
    void find_LU();  
    size_t size;
    Matrix<double> L;
    Matrix<double> U;
    void operator=(const SquareMatrix &second_mat);
    void operator=(const Matrix<double> &second_mat);
};

SquareMatrix::SquareMatrix(const size_t inp_size):  Matrix<double>(inp_size,inp_size), 
                                                    U(inp_size,inp_size), 
                                                    L(inp_size,inp_size), 
                                                    size(inp_size){}

SquareMatrix::SquareMatrix():size(0){}

size_t SquareMatrix::get_size() const{return size;};

double SquareMatrix::get_L_elem(int i, int j) const{return L(i,j);}
double SquareMatrix::get_U_elem(int i, int j) const{return U(i,j);}

void SquareMatrix::find_LU(){
    //Copy matrix in U
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            U(i,j) = elements[i][j];
        }
    }
    for(int k = 0; k < size;++k){
        L(k,k) = 1;
        for(int i = k+1;i<size;++i){
            L(i,k) = U(i,k) / U(k,k);
            for(int j = k; j < size; ++j){
                U(i,j) -= U(k,j)*L(i,k);
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
        ans *= U(i,i);
    }
    return ans;
}

Vector<double> SquareMatrix::find_root(const Vector<double> B){
    if(size == 0){
        std::cerr << "Error: Can't find root with empty matrix" << std::endl;
        throw -1;
    }
    if(size != B.get_size()){
        std::cerr << "Error: Different size of matrix and vector" << std::endl;
        throw -1;
    }

    Vector<double> y(size);
    for(int i = 0;i < size; ++i){
        y(i) = B(i);
        for(int j = 0; j< i;++j){
            y(i) -= y(j)*L(i,j);
        }
        y(i) /= L(i,i);
    }

    Vector<double> x(size);
    for(int i = size - 1;i >= 0; --i){
        x(i) = y(i);
        for(int j = size - 1; j > i;--j){
            x(i) -= x(j)*U(i,j);
        }
        x(i) /= U(i,i);
    }

    return x;
}


SquareMatrix SquareMatrix::get_inverse(){
    SquareMatrix answer(size);
    for(int i = 0; i < size;++i){
        Vector<double> b(size);
        b(i) = 1.0;
        Vector<double> x = find_root(b);
        for(int j = 0; j < size; ++j){
            answer(j,i) = x(j);
        }
    }
    answer.find_LU();
    return answer;
}

void SquareMatrix::operator=(const SquareMatrix &second_mat){
    elements = std::vector<std::vector<double>>(second_mat.get_size(),std::vector<double>(second_mat.get_size()));
    size = second_mat.get_size();

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            elements[i][j] = second_mat(i, j);
        }
    }
}

void SquareMatrix::operator=(const Matrix<double> &second_mat){
    if(second_mat.get_col() != second_mat.get_row()){
        throw std::invalid_argument("different sizes of matrices");
    }

    size = second_mat.get_col();
    elements = std::vector<std::vector<double>>(size,std::vector<double>(size));
    

    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            elements[i][j] = second_mat(i, j);
        }
    }
    
}