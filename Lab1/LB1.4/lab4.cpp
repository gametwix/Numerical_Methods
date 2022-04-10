#include "../matrix.h"
#include <iostream>
#include <cmath>
#include <vector>

Matrix<double> find_U(const Matrix<double> &mat, double eps){
    //find coordinates max abs
    double max_abs = 0;
    int max_row = 0, max_col = 0;
    for(int i = 0; i < mat.get_row();++i){
        for(int j = i + 1; j < mat.get_col();++j){
            if(std::abs(mat(i,j)) > max_abs){
                max_abs = std::abs(mat(i,j));
                max_row = i;
                max_col = j;
            }
        } 
    }

    //find angel phi
    double phi = std::atan((2*mat(max_row,max_col))/(mat(max_row,max_row) - mat(max_col,max_col)))/2;

    //rotation matrix
    Matrix<double> U(mat.get_row(),mat.get_col());

    U(max_row,max_col) = -std::sin(phi);
    U(max_col,max_row) = std::sin(phi);
    U(max_col,max_col) = std::cos(phi);
    U(max_row,max_row) = std::cos(phi);
    for(int i = 0; i < mat.get_row();++i){
        if(i != max_row && i != max_col){
            U(i,i) = 1;
        }
    }

    //matrix for next step
    Matrix<double> new_mat = U.transposition()*mat*U;

    //error search
    double delt = 0.0;
    for(int i = 0; i < mat.get_row();++i){
        for(int j = i + 1; j < mat.get_col();++j){
            delt += new_mat(i,j)*new_mat(i,j);
        } 
    }
    delt = std::sqrt(delt);

    //if error les eps: OK
    if(delt < eps){
        return U;
    }else{
        //find U for next step matrix
        Matrix<double> next_U = find_U(new_mat,eps);
        return U*next_U;
    }
    return new_mat;
}

int main(){
    Matrix<double> mat = {{-7,-5,-9},{-5,5,2},{-9,2,9}};
    double eps = 0.01;
    std::cout << "Matrix:" << std::endl;
    std::cout << mat << std::endl;
    Matrix<double> U = find_U(mat,eps);
    //splitting the matrix into eigenvectors
    std::vector<Vector<double>> vecs(U.get_col(),Vector<double>(U.get_row()));
    std::cout << "Eigenvectors:" << std::endl;
    for(int i = 0; i < U.get_col(); ++i){
        std::cout << "v" << i << std::endl;
        for(int j = 0; j < U.get_row(); ++j){
            vecs[i](j) = U(j,i);
        }
        std::cout << vecs[i] << std::endl;
    }

    //Eigenvalues
    Matrix<double> A = U.transposition()*mat*U;
    std::cout << "Eigenvalues:" ;
    for(int i = 0; i < A.get_col(); ++i){
        std::cout << A(i,i) << " ";
    }

    std::cout << std::endl;
    return 0;
}