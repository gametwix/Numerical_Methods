import numpy as np
from math import sin, cos, exp

class PartialDifferentialEquation:
    def __init__(self, a=1, b=0, c=0,
                 Ut_start_x=lambda t, a, b, c: 0,
                 Ut_finish_x=lambda t, a, b, c: 0,
                 Ux_start_t=lambda x, a, b, c: 0,
                 U_answer=lambda x, t, a, b, c: 0,
                 fi=lambda x, t, a, b, c: 0,
                 ):
        self.a = a
        self.b = b
        self.c = c
        self._f_Ut_start_x = Ut_start_x
        self._f_Ut_finish_x = Ut_finish_x
        self._f_Ux_start_t = Ux_start_t
        self._f_U_answer = U_answer
        self._f_fi = fi

    def Ut_start_x(self, t):
        return self._f_Ut_start_x(t, a=self.a, b=self.b, c=self.c)

    def Ut_finish_x(self, t):
        return self._f_Ut_finish_x(t, a=self.a, b=self.b, c=self.c)

    def Ux_start_t(self, x):
        return self._f_Ux_start_t(x, a=self.a, b=self.b, c=self.c)

    def U_answer(self, x, t):
        return self._f_U_answer(x, t, a=self.a, b=self.b, c=self.c)

    def fi(self, x, t):
        return self._f_fi(x, t, a=self.a, b=self.b, c=self.c)



def Explicit_method(Equation, start_x, finish_x,
                    resolution_x, finish_t,
                    resolution_t
                    ):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    cof_a = Equation.a*tau / h**2
    cof_b = Equation.b*tau / 2*h
    cof_c = Equation.c*tau

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        answer[k+1][0] = Equation.Ut_start_x(k*tau)
        answer[k+1][resolution_x - 1] = Equation.Ut_finish_x(k*tau)
        for i in range(1, resolution_x - 1):
            answer[k+1][i] += (cof_c + 1 - 2*cof_a) * answer[k][i]
            answer[k+1][i] += (cof_a - cof_b) * answer[k][i-1]
            answer[k+1][i] += (cof_a + cof_b) * answer[k][i+1]
            answer[k+1][i] += Equation.fi(i*h + start_x, tau*k)
    return answer


def Implicit_method(Equation, start_x, finish_x,
                    resolution_x, finish_t,
                    resolution_t
                    ):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    cof_a = Equation.a*tau / h**2
    cof_b = Equation.b*tau / 2*h
    cof_c = Equation.c*tau

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        answer[k+1][0] = Equation.Ut_start_x(k*tau)
        answer[k+1][resolution_x - 1] = Equation.Ut_finish_x(k*tau)
        new_line_mat = np.zeros((resolution_x-2, resolution_x-2))
        new_line_vec = np.zeros(resolution_x-2)
        for i in range(resolution_x-2):
            new_line_mat[i][i] = (2*cof_a + cof_c + 1)
            new_line_vec[i] = answer[k][i+1] + \
                Equation.fi((i+1) * h + start_x, tau * (k+1))

            if i == 0:
                new_line_mat[i][i+1] = -(cof_b + cof_a)
                new_line_vec[i] += (cof_a - cof_b)*answer[k+1][0]
            elif i == resolution_x-3:
                new_line_mat[i][i-1] = cof_b - cof_a
                new_line_vec[i] += (cof_a + cof_b)*answer[k+1][resolution_x - 1]
            else:
                new_line_mat[i][i+1] = -(cof_b + cof_a)
                new_line_mat[i][i-1] = cof_b - cof_a

        line = np.linalg.solve(new_line_mat, new_line_vec)
        answer[k+1][1:-1] = line
    return answer


def Explicit_Implicit_method(Equation, start_x, finish_x,
                             resolution_x, finish_t,
                             resolution_t, teta=0.5
                             ):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    cof_a_imp = Equation.a*tau*teta / h**2
    cof_b_imp = Equation.b*tau*teta / 2*h
    cof_c_imp = Equation.c*tau*teta
    cof_a_exp = Equation.a*tau*(1-teta) / h**2
    cof_b_exp = Equation.b*tau*(1-teta) / 2*h
    cof_c_exp = Equation.c*tau*(1-teta)

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        answer[k+1][0] = Equation.Ut_start_x(k*tau)
        answer[k+1][resolution_x - 1] = Equation.Ut_finish_x(k*tau)
        new_line_mat = np.zeros((resolution_x-2, resolution_x-2))
        new_line_vec = np.zeros(resolution_x-2)
        for i in range(resolution_x-2):
            new_line_mat[i][i] = (2*cof_a_imp + cof_c_imp + 1)
            new_line_vec[i] = teta*Equation.fi((i+1) * h + start_x, tau * (k+1))\
                + (cof_c_exp + 1 - 2*cof_a_exp) * answer[k][i+1]\
                + (cof_a_exp - cof_b_exp) * answer[k][i]\
                + (cof_a_exp + cof_b_exp) * answer[k][i+2]\
                + (1 - teta)*Equation.fi((i+1) * h + start_x, tau * k)

            if i == 0:
                new_line_mat[i][i+1] = -(cof_b_imp+cof_a_imp)
                new_line_vec[i] += (cof_a_imp-cof_b_imp) * answer[k+1][0]
            elif i == resolution_x-3:
                new_line_mat[i][i-1] = cof_b_imp - cof_a_imp
                new_line_vec[i] += (cof_a_imp+cof_b_imp) * answer[k+1][resolution_x - 1]
            else:
                new_line_mat[i][i+1] = -(cof_b_imp + cof_a_imp)
                new_line_mat[i][i-1] = cof_b_imp - cof_a_imp

        line = np.linalg.solve(new_line_mat, new_line_vec)
        answer[k+1][1:-1] = line
    return answer


    