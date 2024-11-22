extern "C" {
double* Volterra_int_eq_sec_kind_solv(int N, double a, double b,
                                            double (*K)(double x, double y),
                                            double (*f)(double x),
                                            double beta);
}
