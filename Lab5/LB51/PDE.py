import numpy as np

def trid_alg(matrix,vec):
    s = len(vec)
    ans = np.zeros(s)
    p = np.zeros(s)
    q = np.zeros(s)
    
    p[0] = -matrix[0][1] / matrix[0][0]
    q[0] = vec[0] / matrix[0][0]
    
    for i in range(1,s - 1):
        p[i] = -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i-1]*p[i-1])
        q[i] = (vec[i] - matrix[i][i-1]*q[i-1]) / (matrix[i][i] + matrix[i][i-1]*p[i-1])
        
    p[s-1] = 0
    q[s-1] = (vec[s-1] - matrix[s-1][s-2]*q[s-2]) / (matrix[s-1][s-1] + matrix[s-1][s-2]*p[s-2])
    
    ans[s-1] = q[s-1]
    
    for i in range(s - 2, -1, -1):
        ans[i] = p[i]*ans[i+1] + q[i]
        
    return ans


class PartialDifferentialEquation:
    def __init__(
        self,
        a=1,
        b=0,
        c=0,
        Ut_start_x=lambda t, a, b, c: 0,
        Ut_finish_x=lambda t, a, b, c: 0,
        Ux_start_t=lambda x, a, b, c: 0,
        U_answer=lambda x, t, a, b, c: 0,
        fi=lambda x, t, a, b, c: 0,
        alf0=1,
        bet0=0,
        alfn=1,
        betn=0,
    ):
        self.a = a
        self.b = b
        self.c = c
        self.alf0 = alf0
        self.bet0 = bet0
        self.alfn = alfn
        self.betn = betn
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


def Explicit_method(Equation,
                    start_x,
                    finish_x,
                    resolution_x,
                    finish_t,
                    resolution_t,
                    method=3
                    ):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h * i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    cof_a = Equation.a * tau / h**2
    cof_b = Equation.b * tau / 2 * h
    cof_c = Equation.c * tau
    #print(Equation.a*tau/(h*h))

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        for i in range(1, resolution_x - 1):
            answer[k + 1][i] += ((cof_c + 1 - 2 * cof_a) * answer[k][i]
                                 + (cof_a - cof_b) * answer[k][i - 1]
                                 + (cof_a + cof_b) * answer[k][i + 1]
                                 + Equation.fi(i * h + start_x, tau * (k+1)))
            
        if method == 1:
            b0 = Equation.bet0 - Equation.alf0/h
            c0 = Equation.alf0 / h
            d0 = Equation.Ut_start_x((k+1) * tau)
            an = -Equation.alfn/h
            bn = Equation.betn + Equation.alfn/h
            dn = Equation.Ut_finish_x((k+1) * tau)
            answer[k + 1][0] = (d0 - answer[k+1][1]*c0) / b0
            answer[k + 1][resolution_x - 1] = (dn - answer[k + 1][resolution_x - 2]*an) / bn
        elif method == 2:
            b0 = Equation.bet0 - (3*Equation.alf0) / (2*h)
            c0 = 2*Equation.alf0 / h
            e0 = - Equation.alf0 / (2*h)
            bn = Equation.betn + (3*Equation.alfn) / (2*h)
            cn = - 2*Equation.alfn / h
            en = + Equation.alfn / (2*h)
            answer[k + 1][0] = (Equation.Ut_start_x((k+1) * tau)
                                - c0*answer[k+1][1]
                                - e0*answer[k+1][2]
                                ) / b0
            answer[k + 1][resolution_x - 1] = (Equation.Ut_finish_x((k+1) * tau)
                                               - cn*answer[k+1][resolution_x - 2]
                                               - en*answer[k+1][resolution_x - 3]) / bn
        elif method == 3:
            answer[k + 1][0] = (h*Equation.alf0*answer[k][0]/tau
                                - Equation.Ut_start_x((k+1) * tau)
                                * (2*Equation.a - Equation.b*h)
                                - Equation.alf0*h*Equation.fi(0, (k+1) * tau))
            answer[k + 1][0] -= -2*Equation.a*answer[k+1][1]*Equation.alf0/h
            answer[k + 1][0] /= (2*Equation.a*Equation.alf0 / h
                                    + h*Equation.alf0 / tau
                                    - Equation.c * h*Equation.alf0
                                    - Equation.bet0*(2*Equation.a - Equation.b*h)
                                )
            
            answer[k + 1][resolution_x - 1] = (h*Equation.alfn*answer[k][resolution_x - 1]/tau
                                                + Equation.Ut_finish_x((k+1) * tau)
                                                * (2*Equation.a - Equation.b*h)
                                                + Equation.alfn*h*Equation.fi(finish_x, (k+1) * tau))
            answer[k + 1][resolution_x - 1] -= -2*Equation.a*Equation.alfn*answer[k+1][resolution_x - 2]/h
            answer[k + 1][resolution_x - 1] /= (2*Equation.a*Equation.alfn / h
                    + h*Equation.alfn / tau
                    - Equation.c * h*Equation.alfn
                    + Equation.betn*(2*Equation.a - Equation.b*h))
    return answer


def Implicit_method(Equation,
                    start_x,
                    finish_x,
                    resolution_x,
                    finish_t,
                    resolution_t,
                    method=3
                    ):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h * i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    cof_a = Equation.a * tau / h**2
    cof_b = Equation.b * tau / 2 * h
    cof_c = Equation.c * tau

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        new_line_mat = np.zeros((resolution_x, resolution_x))
        new_line_vec = np.zeros(resolution_x)
        if method == 1:
            b0 = Equation.bet0 - Equation.alf0/h
            c0 = Equation.alf0 / h
            d0 = Equation.Ut_start_x((k+1) * tau) 
            an = -Equation.alfn/h
            bn = Equation.betn + Equation.alfn/h
            dn = Equation.Ut_finish_x((k+1) * tau)
            new_line_mat[0][0] = b0
            new_line_mat[0][1] = c0
            new_line_mat[resolution_x - 1][resolution_x - 1] = bn
            new_line_mat[resolution_x - 1][resolution_x - 2] = an
            new_line_vec[0] = d0
            new_line_vec[resolution_x - 1] = dn
        elif method == 2:
            b0 = Equation.bet0 - (3*Equation.alf0) / (2*h)
            c0 = 2*Equation.alf0 / h
            e0 = - Equation.alf0 / (2*h)
            bn = Equation.betn + (3*Equation.alfn) / (2*h)
            cn = - 2*Equation.alfn / h
            en = Equation.alfn / (2*h)
            new_line_mat[0][0] = b0
            new_line_mat[0][1] = c0
            new_line_mat[0][2] = e0
            new_line_mat[resolution_x - 1][resolution_x - 1] = bn
            new_line_mat[resolution_x - 1][resolution_x - 2] = cn
            new_line_mat[resolution_x - 1][resolution_x - 3] = en
            new_line_vec[0] = Equation.Ut_start_x((k+1) * tau)
            new_line_vec[resolution_x - 1] = Equation.Ut_finish_x((k+1) * tau)
        elif method == 3:
            new_line_mat[0][0] = (2*Equation.a*Equation.alf0 / h
                                + h*Equation.alf0 / tau
                                - Equation.c * h*Equation.alf0
                                - Equation.bet0*(2*Equation.a - Equation.b*h)
                                )
            new_line_mat[0][1] = -2*Equation.a*Equation.alf0/h
            
            new_line_mat[resolution_x - 1][resolution_x - 1] = \
                (2*Equation.a*Equation.alfn / h
                + h*Equation.alfn / tau
                - Equation.c * h*Equation.alfn
                + Equation.betn*(2*Equation.a - Equation.b*h))
            new_line_mat[resolution_x - 1][resolution_x - 2] = -2*Equation.a*Equation.alfn/h
            
            new_line_vec[0] = (h*Equation.alf0*answer[k+1][0]/tau
                            - Equation.Ut_start_x((k+1) * tau)
                            * (2*Equation.a - Equation.b*h)
                            - Equation.alf0*h*Equation.fi(0, (k+1) * tau))
            new_line_vec[resolution_x - 1] = (h*Equation.alfn*answer[k+1][resolution_x - 1]/tau
                                            + Equation.Ut_finish_x((k+1) * tau)
                                            * (2*Equation.a - Equation.b*h)
                                            + Equation.alfn*h*Equation.fi(finish_x, (k+1) * tau))
        
        for i in range(1, resolution_x - 1):
            new_line_mat[i][i - 1] = cof_b - cof_a
            new_line_mat[i][i] = 2 * cof_a + cof_c + 1
            new_line_mat[i][i + 1] = -(cof_b + cof_a)

            new_line_vec[i] = answer[k][i] + Equation.fi(
                i * h + start_x, tau * (k + 1)
            )

        if method == 2:
            k0 = new_line_mat[0][2] / new_line_mat[1][2]
            kn = (new_line_mat[resolution_x - 1][resolution_x - 3]
                  / new_line_mat[resolution_x - 2][resolution_x - 3])
            for i in range(resolution_x):
                new_line_mat[0][i] -= new_line_mat[1][i] * k0
                new_line_mat[resolution_x - 1][i] -= new_line_mat[resolution_x - 2][i] * kn
            new_line_vec[0] -= new_line_vec[1] * k0
            new_line_vec[resolution_x - 1] -= new_line_vec[resolution_x - 2] * kn

        line = trid_alg(new_line_mat, new_line_vec)
        answer[k + 1] = line
    return answer


def Explicit_Implicit_method(Equation,
                             start_x,
                             finish_x,
                             resolution_x,
                             finish_t,
                             resolution_t,
                             teta=0.5,
                             method=3):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h * i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))

    #  From Implicit
    cof_a_imp = Equation.a * tau * teta / h**2
    cof_b_imp = Equation.b * tau * teta / 2 * h
    cof_c_imp = Equation.c * tau * teta
    #  From Explicit
    cof_a_exp = Equation.a * tau * (1 - teta) / h**2
    cof_b_exp = Equation.b * tau * (1 - teta) / 2 * h
    cof_c_exp = Equation.c * tau * (1 - teta)

    answer[0, :] = [Equation.Ux_start_t(elem) for elem in x]
    for k in range(resolution_t - 1):
        #  From Implicit
        new_line_mat = np.zeros((resolution_x, resolution_x))
        new_line_vec = np.zeros(resolution_x)
        if method == 1:
            b0 = Equation.bet0 - Equation.alf0/h
            c0 = Equation.alf0 / h
            d0 = Equation.Ut_start_x((k+1) * tau) 
            an = -Equation.alfn/h
            bn = Equation.betn + Equation.alfn/h
            dn = Equation.Ut_finish_x((k+1) * tau)
            new_line_mat[0][0] = b0
            new_line_mat[0][1] = c0
            new_line_mat[resolution_x - 1][resolution_x - 1] = bn
            new_line_mat[resolution_x - 1][resolution_x - 2] = an
            new_line_vec[0] = d0
            new_line_vec[resolution_x - 1] = dn
        elif method == 2:
            b0 = Equation.bet0 - (3*Equation.alf0) / (2*h)
            c0 = 2*Equation.alf0 / h
            e0 = - Equation.alf0 / (2*h)
            bn = Equation.betn + (3*Equation.alfn) / (2*h)
            cn = - 2*Equation.alfn / h
            en =  Equation.alfn / (2*h)
            new_line_mat[0][0] = b0
            new_line_mat[0][1] = c0
            new_line_mat[0][2] = e0
            new_line_mat[resolution_x - 1][resolution_x - 1] = bn
            new_line_mat[resolution_x - 1][resolution_x - 2] = cn
            new_line_mat[resolution_x - 1][resolution_x - 3] = en
            new_line_vec[0] = Equation.Ut_start_x((k+1) * tau)
            new_line_vec[resolution_x - 1] = Equation.Ut_finish_x((k+1) * tau)
        elif method == 3:
            new_line_mat[0][0] = (2*Equation.a*Equation.alf0 / h
                                + h*Equation.alf0 / tau
                                - Equation.c * h*Equation.alf0
                                - Equation.bet0*(2*Equation.a - Equation.b*h)
                                )
            new_line_mat[0][1] = -2*Equation.a*Equation.alf0/h
            
            new_line_mat[resolution_x - 1][resolution_x - 1] = \
                (2*Equation.a*Equation.alfn / h
                + h*Equation.alfn / tau
                - Equation.c * h*Equation.alfn
                + Equation.betn*(2*Equation.a - Equation.b*h))
            new_line_mat[resolution_x - 1][resolution_x - 2] = -2*Equation.a*Equation.alfn/h
            
            new_line_vec[0] = (h*Equation.alf0*answer[k+1][0]/tau
                            - Equation.Ut_start_x((k+1) * tau)
                            * (2*Equation.a - Equation.b*h)
                            - Equation.alf0*h*Equation.fi(0, (k+1) * tau))
            new_line_vec[resolution_x - 1] = (h*Equation.alfn*answer[k+1][resolution_x - 1]/tau
                                            + Equation.Ut_finish_x((k+1) * tau)
                                            * (2*Equation.a - Equation.b*h)
                                            + Equation.alfn*h*Equation.fi(finish_x, (k+1) * tau))
        for i in range(1, resolution_x - 1):
            #  From Implicit
            new_line_mat[i][i - 1] = cof_b_imp - cof_a_imp
            new_line_mat[i][i] = 2 * cof_a_imp + cof_c_imp + 1
            new_line_mat[i][i + 1] = -(cof_b_imp + cof_a_imp)
            
            new_line_vec[i] = (
                #  From Implicit
                teta * Equation.fi(i * h + start_x, tau * (k + 1))
                #  From Explicit
                + (cof_c_exp + 1 - 2 * cof_a_exp) * answer[k][i]
                + (cof_a_exp - cof_b_exp) * answer[k][i - 1]
                + (cof_a_exp + cof_b_exp) * answer[k][i + 1]
                + (1 - teta) * Equation.fi(i * h + start_x, tau * k)
            )
            
        if method == 2:
            k0 = new_line_mat[0][2] / new_line_mat[1][2] 
            kn = (new_line_mat[resolution_x - 1][resolution_x - 3]
                  / new_line_mat[resolution_x - 2][resolution_x - 3])
            for i in range(resolution_x):
                new_line_mat[0][i] -= new_line_mat[1][i] * k0
                new_line_mat[resolution_x - 1][i] -= new_line_mat[resolution_x - 2][i] * kn
            new_line_vec[0] -= new_line_vec[1] * k0
            new_line_vec[resolution_x - 1] -= new_line_vec[resolution_x - 2] * kn
            
        line = trid_alg(new_line_mat, new_line_vec)
        answer[k + 1] = line
    return answer
