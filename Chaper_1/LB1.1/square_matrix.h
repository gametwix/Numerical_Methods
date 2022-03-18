#pragma once
#include <iostream>
#include <vector>

class SquareMatrix{
private:
    int size;
    std::vector<std::vector<double>> elements;
    std::vector<std::vector<double>> L;
    std::vector<std::vector<double>> U;

    void find_LU();
public:
    SquareMatrix();
    SquareMatrix(const int _size);
    SquareMatrix(const std::vector<std::vector<double>> _elements);
    SquareMatrix(const SquareMatrix &_matrix);

    int get_size();
    double get_elem(int i, int j);
    double get_L_elem(int i, int j);
    double get_U_elem(int i, int j);
    double get_det();
    SquareMatrix get_inverse();

    void set_elem(int i,int j,double elem);


    std::vector<double> find_root(const std::vector<double> B);
};