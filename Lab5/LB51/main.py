import numpy as np
import matplotlib.pyplot as plt


def Ut_start_x(t, a = 1):
    return 0

def Ut_finish_x(t, a = 1):
    return 0

def Ux_start_t(x, a = 1):
    return np.sin(x)

def U_answer(x,t, a = 1):
    return np.exp(-a*t)*np.sin(x)

def Explicit_method(start_x, finish_x, resolution_x, finish_t, resolution_t):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))
    # Избавится потом об этого переменными окружения, параметрами функции или завертыванием в класс
    a = 1
    b = 0
    c= 0
    fi = lambda x, t: 0
    #----
    
    cof_a = a*tau / h**2
    cof_b = b*tau / 2*h
    cof_c = c*tau
    
    answer[0,:] = Ux_start_t(x, a)
    for k in range(resolution_t - 1):
        answer[k+1][0] = Ut_start_x(k*tau,a)
        answer[k+1][resolution_x - 1] = Ut_finish_x(k*tau,a)
        for i in range(1,resolution_x - 1):            
            answer[k+1][i] += (cof_c + 1 - 2*cof_a) * answer[k][i]
            answer[k+1][i] += (cof_a - cof_b) * answer[k][i-1]
            answer[k+1][i] += (cof_a + cof_b) * answer[k][i+1]
            answer[k+1][i] += fi(i * h + start_x, tau * k)
    return answer

def Implicit_method(start_x, finish_x, resolution_x, finish_t, resolution_t):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))
    # Избавится потом об этого переменными окружения, параметрами функции или завертыванием в класс
    a = 1
    b = 0
    c= 0
    fi = lambda x, t: 0
    #----
    
    cof_a = a*tau / h**2
    cof_b = b*tau / 2*h
    cof_c = c*tau
    
    answer[0,:] = Ux_start_t(x, a)
    for k in range(resolution_t - 1):
        answer[k+1][0] = Ut_start_x(k*tau,a)
        answer[k+1][resolution_x - 1] = Ut_finish_x(k*tau,a)
        new_line_mat = np.zeros((resolution_x-2, resolution_x-2))
        new_line_vec = np.zeros(resolution_x-2)
        for i in range(resolution_x-2):
            new_line_mat[i][i] = (2*cof_a + cof_c + 1)
            new_line_vec[i] = answer[k][i+1] + fi((i+1) * h + start_x, tau * (k+1))
            
            if i == 0:
                new_line_mat[i][i+1] = -(cof_b + cof_a)
                new_line_vec[i] += (cof_a - cof_b)*answer[k+1][0]
            elif i == resolution_x-3:
                new_line_mat[i][i-1] = cof_b - cof_a
                new_line_vec[i] += (cof_a + cof_b)*answer[k+1][resolution_x - 1]
            else:
                new_line_mat[i][i+1] = -(cof_b + cof_a)
                new_line_mat[i][i-1] = cof_b - cof_a
                
        line = np.linalg.solve(new_line_mat,new_line_vec)
        answer[k+1][1:-1] = line
    return answer

def Explicit_Implicit_method(start_x, finish_x, resolution_x, finish_t, resolution_t, teta = 0.5):
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = finish_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    answer = np.zeros((resolution_t, resolution_x))
    # Избавится потом об этого переменными окружения, параметрами функции или завертыванием в класс
    a = 1
    b = 0
    c= 0
    fi = lambda x, t: 0
    #----
    
    cof_a_imp = a*tau*teta / h**2
    cof_b_imp = b*tau*teta / 2*h
    cof_c_imp = c*tau*teta
    cof_a_exp = a*tau*(1-teta) / h**2
    cof_b_exp = b*tau*(1-teta) / 2*h
    cof_c_exp = c*tau*(1-teta)
    
    answer[0,:] = Ux_start_t(x, a)
    for k in range(resolution_t - 1):
        answer[k+1][0] = Ut_start_x(k*tau,a)
        answer[k+1][resolution_x - 1] = Ut_finish_x(k*tau,a)
        new_line_mat = np.zeros((resolution_x-2, resolution_x-2))
        new_line_vec = np.zeros(resolution_x-2)
        for i in range(resolution_x-2):
            new_line_mat[i][i] = (2*cof_a_imp + cof_c_imp + 1)
            new_line_vec[i] = teta*fi((i+1) * h + start_x, tau * (k+1))\
                + (cof_c_exp + 1 - 2*cof_a_exp) * answer[k][i+1]\
                + (cof_a_exp - cof_b_exp) * answer[k][i]\
                + (cof_a_exp + cof_b_exp) * answer[k][i+2]\
                + (1 - teta)*fi((i+1) * h + start_x, tau * k)

            if i == 0:
                new_line_mat[i][i+1] = -(cof_b_imp + cof_a_imp)
                new_line_vec[i] += (cof_a_imp - cof_b_imp)*answer[k+1][0]
            elif i == resolution_x-3:
                new_line_mat[i][i-1] = cof_b_imp - cof_a_imp
                new_line_vec[i] += (cof_a_imp + cof_b_imp)*answer[k+1][resolution_x - 1]
            else:
                new_line_mat[i][i+1] = -(cof_b_imp + cof_a_imp)
                new_line_mat[i][i-1] = cof_b_imp - cof_a_imp
                
        line = np.linalg.solve(new_line_mat,new_line_vec)
        answer[k+1][1:-1] = line
    return answer
                    
    
if __name__ == '__main__':
    start_x = 0
    finish_x = np.pi
    resolution_x = 30
    h = (finish_x - start_x)  / (resolution_x - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    y_ans = U_answer(x,0.5)
    plt.plot(x,y_ans)
    explicit = Explicit_Implicit_method(start_x, finish_x, resolution_x, 1, 11)
    
    plt.plot(x,explicit[5])
    plt.show()
    