#include "include/matrix.h"
#include "include/square_matrix.h"
#include "include/volterra_solv.h"
extern "C" {
double* Volterra_int_eq_sec_kind_solv(int N, double a, double b,
                                            double (*K)(double x, double y),
                                            double (*f)(double x),
                                            double beta){
    //шаг интегрирования методом трапеции (шаг между узлами в квадратурной сумме)
    double h = (b - a) / (N - 1);
    //Веса квадратурной суммы (интегралл методом трапеции)
    Vector<double> weight(N);
    for(int i = 0; i < N; ++i){
        if(i == 0 || i == (N-1)){
            weight(i) = h/2;
        }else{
            weight(i) = h;
        }
    }

    //Коэфициэнты элеметов под интеграллом для создания системы уравнений
    SquareMatrix C(N);
    for(int i = 0; i < N; ++i){
        //Для заполнения не нижнетреугольной матрицы, а полной нужно заменить j<=i на j < N
        //это позволит решать уравнения Фредгольма 
        for(int j = 0; j <= i; ++j){
            C(i,j) = beta*weight(j)*K(a+i*h,a+j*h);
        }
    }

    //Единичная матрица
    SquareMatrix I(N);
    for(int i = 0; i < N; ++i){
        I(i,i) = 1.;
    }

    //Коэфиценты левой части системы уравнений
    SquareMatrix B(N);
    B = I - C;

    //Коэфиценты правой части системы уравнений
    Vector<double> R(N);
    for(int i = 0; i < N; ++i){
        //получаем значения функции с шагом h
        R(i) = f(a+i*h);
    }

    B.find_LU();
    //Решаем систему уравнений (В моем случаем методом Гауса)
    Vector<double> ans = B.find_root(R);

    double* ans_m = new double[N];
    for(int i = 0; i < N; ++i){
        ans_m[i] = ans(i);
    }
    return ans_m;
}
}