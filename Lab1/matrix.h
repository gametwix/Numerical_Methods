#pragma once
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

template <class T>
class Matrix{
public:
    Matrix():columns(0),rows(0){};
    Matrix(const size_t row_size,const size_t col_size);
    Matrix(std::initializer_list<std::initializer_list<double>> elems);

    T operator()(const size_t row,const size_t col) const;
    T& operator()(const size_t row,const size_t col);
    Matrix operator+(const Matrix<T> &second_mat);
    Matrix operator-(const Matrix<T> &second_mat);
    Matrix operator*(const Matrix<T> &second_mat);
    void operator=(const Matrix<T> &second_mat);
    size_t get_col() const {return columns;}
    size_t get_row() const {return rows;}

    friend std::ostream& operator<<(std::ostream &os, const Matrix<T> &matrix){
        os << std::fixed << std::setprecision(3);
        for(size_t i = 0; i < matrix.rows;++i){
            for(size_t j = 0; j < matrix.columns;++j){
                os << std::setw(7) << matrix(i,j);
            }
            os << std::endl;
        }
        os << std::defaultfloat << std::setprecision(6);
        return os;
    }

    std::vector<std::vector<T>> elements;
    size_t columns, rows;
};

template <class T>
Matrix<T>::Matrix(size_t row_size,size_t col_size):columns(col_size),rows(row_size){
    elements.resize(rows);
    for(int i = 0; i < rows;++i){
        elements[i].resize(columns);
    }
}

template <class T>
T& Matrix<T>::operator()(const size_t row,const size_t col){
    if(col >= columns || row >= rows){
        throw std::out_of_range("arguments out of range");
    }
    return elements[row][col];
}

template <class T>
T Matrix<T>::operator()(const size_t row,const size_t col) const{
    if(col >= columns || row >= rows){
        throw std::out_of_range("arguments out of range");
    }
    return elements[row][col];
}

template <class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &second_mat){
    if(columns != second_mat.get_col() || rows != second_mat.get_row()){
        throw std::invalid_argument("different sizes of matrices");
    }

    Matrix<T> answer(rows,columns);

    for(int i = 0;i < rows; ++i){
        for(int j = 0; j < columns; ++j){
            answer(i,j) = elements[i][j] + second_mat(i,j);
        }
    }

    return answer;
}

template <class T>
void Matrix<T>::operator=(const Matrix<T> &second_mat){
    elements = std::vector<std::vector<T>>(second_mat.get_row(),std::vector<T>(second_mat.get_col()));
    rows = second_mat.get_row();
    columns = second_mat.get_col();
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < columns; ++j){
            elements[i][j] = second_mat(i, j);
        }
    }
}

template <class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &second_mat){
    if(columns != second_mat.get_col() || rows != second_mat.get_row()){
        throw std::invalid_argument("different sizes of matrices");
    }

    Matrix<T> answer(rows,columns);

    for(size_t i = 0;i < rows; ++i){
        for(size_t j = 0; j < columns; ++j){
            answer(i,j) = elements[i][j] - second_mat(i,j);
        }
    }

    return answer;
}

template <class T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &second_mat){
    if(columns != second_mat.get_row()){
        throw std::invalid_argument("different sizes of matrices");
    }

    Matrix<T> answer(rows, second_mat.get_col());
    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < second_mat.get_col(); ++j){
            answer(i,j) = 0;
        }
    }

    for(size_t i = 0; i < rows; ++i){
        for(size_t j = 0; j < second_mat.get_col(); ++j){
            for(size_t k = 0; k < columns; ++k){
                answer(i,j) += elements[i][k]*second_mat(k,j);
            }
        }
    }
    return answer;
}

template <class T>
Matrix<T>::Matrix(std::initializer_list<std::initializer_list<double>> elems){
    elements.resize(elems.size());
    size_t i = 0;
    size_t j = 0;
    for(auto list:elems){
        j = 0;
        elements[i].resize(list.size());
        for(auto elem: list){
            elements[i][j] = elem;
            j++;
        }
        i++;
    }
    rows = i;
    columns = j;
}



template <class T>
class Vector: public Matrix<T>{
    using Matrix<T>::rows;
    using Matrix<T>::columns;
    using Matrix<T>::elements;
public:
    Vector(const size_t inp_size): size(inp_size),Matrix<T>::Matrix(inp_size,1){}
    Vector(std::initializer_list<T> list):Vector(list.size()){
        size_t i = 0;
        for(auto elem:list){
            elements[i][0] = elem;
            i++;
        }
    }
    Vector():size(0){}

    size_t get_size() const{return size;}

    T& operator()(const size_t num){
        return Matrix<T>::operator()(num,0);
    }

    T operator()(const size_t num)const{
        return Matrix<T>::operator()(num,0);
    }

    void operator=(const Matrix<T> &second_mat){
        if(second_mat.get_col() > 1){
            throw std::invalid_argument("different sizes of matrices");
        }
        elements = std::vector<std::vector<T>>(second_mat.get_row(),std::vector<T>(1));
        rows = second_mat.get_row();
        size = second_mat.get_row();
        columns = second_mat.get_col();
        for(int i = 0; i < rows; ++i){
            elements[i][0] = second_mat(i, 0);
        }
    }
private:
    size_t size;
};
