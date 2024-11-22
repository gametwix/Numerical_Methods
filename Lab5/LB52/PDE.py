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
        b=1,
        c=-1,
        e=3,
        Ut_start_x=lambda t, a, b, c: 0,
        Ut_finish_x=lambda t, a, b, c: 0,
        Ux_start_t=lambda x, a, b, c: 0,
        Ux_start_t_dt=lambda x, a, b, c: 0,
        U_answer=lambda x, t, a, b, c: 0,
        Ux_start_t_dx=lambda x, a, b, c: 0,
        Ux_start_t_dxdx=lambda x, a, b, c: 0,
        fi=lambda x, t, a, b, c: 0,
        alf0=1,
        bet0=-2,
        alfn=1,
        betn=-2,
    ):
        self.a = a
        self.b = b
        self.c = c
        self.e = e
        
        self.alf0 = alf0
        self.bet0 = bet0
        self.alfn = alfn
        self.betn = betn 
        
        # U(0, t) and U(l, t)
        self._f_Ut_start_x = Ut_start_x
        self._f_Ut_finish_x = Ut_finish_x
        
        # U(x, 0) and dU(x, 0) / dx
        self._f_Ux_start_t = Ux_start_t
        self._f_Ux_start_t_dt = Ux_start_t_dt
        self._f_Ux_start_t_dx = Ux_start_t_dx
        self._f_Ux_start_t_dxdx = Ux_start_t_dxdx
        
        # f(x, t)
        self._f_fi = fi
        self._f_U_answer = U_answer
        
    
    # U(0, t) and U(l, t)
    def Ut_start_x(self, t):
        return self._f_Ut_start_x(t, a=self.a, b=self.b, c=self.c)

    def Ut_finish_x(self, t):
        return self._f_Ut_finish_x(t, a=self.a, b=self.b, c=self.c)
    
    # U(x, 0) and dU(x, 0) / dx
    def Ux_start_t(self, x):
        return self._f_Ux_start_t(x, a=self.a, b=self.b, c=self.c)
    
    def Ux_start_t_dt(self, x):
        return self._f_Ux_start_t_dt(x, a=self.a, b=self.b, c=self.c)
    
    def Ux_start_t_dx(self, x):
        return self._f_Ux_start_t_dx(x, a=self.a, b=self.b, c=self.c)
    
    def Ux_start_t_dxdx(self, x):
        return self._f_Ux_start_t_dxdx(x, a=self.a, b=self.b, c=self.c)

    def U_answer(self, x, t):
        return self._f_U_answer(x, t, a=self.a, b=self.b, c=self.c)

    def fi(self, x, t):
        return self._f_fi(x, t, a=self.a, b=self.b, c=self.c)



        
    
        

