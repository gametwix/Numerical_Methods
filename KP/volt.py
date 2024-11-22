import ctypes
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import parser
from math import sin, cos, exp


lib = ctypes.cdll.LoadLibrary('libVolt.so')


def K(x,y):
    global K_str
    formula = K_str
    code = parser.expr(formula).compile()
    return eval(code)

def f(x):
    global f_str
    formula = f_str
    code = parser.expr(formula).compile()
    return eval(code)

def real_fn(x):
    global ans_str
    formula = ans_str
    code = parser.expr(formula).compile()
    return eval(code)

def Volterra(N,a,b,K,f,beta):
    lib.Volterra_int_eq_sec_kind_solv.restype = ctypes.POINTER(ctypes.c_double * N)
    ft_2 = ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double,ctypes.c_double)
    ft_1 = ctypes.CFUNCTYPE(ctypes.c_double,ctypes.c_double)
    lib.Volterra_int_eq_sec_kind_solv.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.c_double,ft_2,ft_1,ctypes.c_double]
    Array = lib.Volterra_int_eq_sec_kind_solv(ctypes.c_int(N),ctypes.c_double(a),ctypes.c_double(b),ft_2(K),ft_1(f),ctypes.c_double(beta))
    return list(Array.contents)


N = 50
a = 0
b = 2
beta = 1
K_str = "(1 + x**2)/(1+y**2)"
f_str = "1 + x**2"
ans_str = "(1+x**2)*exp(x)"

ans = Volterra(N,a,b,K,f,beta)

x_ans = np.linspace(a, b, N)
y_ans = [real_fn(xi) for xi in x_ans]

fig = plt.figure()
fig.subplots_adjust(right=0.7)
ax = fig.add_subplot()
ax.grid()
real, = ax.plot(x_ans, y_ans)
my_plot, = ax.plot(x_ans, ans)


def update_function(param_name, param_val):
    global N,a,b,beta,K_str,f_str,ans_str
    if param_name == 'N':
        N = param_val
    elif param_name == 'a':
        a = param_val
    elif param_name == 'b':
        b = param_val
    elif param_name == 'beta':
        beta = param_val
    elif param_name == 'K_str':
        K_str = param_val
    elif param_name == 'f_str':
        f_str = param_val
    elif param_name == 'ans':
        ans_str = param_val
    x_ans = np.linspace(a, b, N)
    y_ans = [real_fn(xi) for xi in x_ans]
    ans = Volterra(N,a,b,K,f,beta)
    real.set_xdata(x_ans)
    real.set_ydata(y_ans)
    my_plot.set_xdata(x_ans)
    my_plot.set_ydata(ans)
    ax.relim()
    ax.autoscale_view()

def submit_fn_a(value):
    update_function("a",float(value))
    plt.draw()
axbox_a = fig.add_axes([0.75, 0.93, 0.1, 0.05])
text_box_a = TextBox(axbox_a, "a ")
text_box_a.on_submit(submit_fn_a)
text_box_a.set_val(a)

def submit_fn_b(value):
    update_function("b",float(value))
    plt.draw()
axbox_b = fig.add_axes([0.89, 0.93, 0.1, 0.05])
text_box_b = TextBox(axbox_b, "b ")
text_box_b.on_submit(submit_fn_b)
text_box_b.set_val(b)

def submit_fn_N(value):
    update_function("N",int(value))
    plt.draw()
axbox_N = fig.add_axes([0.75, 0.86, 0.1, 0.05])
text_box_N = TextBox(axbox_N, "N ")
text_box_N.on_submit(submit_fn_N)
text_box_N.set_val(N)

def submit_fn_beta(value):
    update_function("beta",float(value))
    plt.draw()
axbox_beta = fig.add_axes([0.89, 0.86, 0.1, 0.05])
text_box_beta = TextBox(axbox_beta, "B ")
text_box_beta.on_submit(submit_fn_beta)
text_box_beta.set_val(beta)

def submit_fn_K(value):
    update_function("K_str",str(value))
    plt.draw()
axbox_K = fig.add_axes([0.79, 0.79, 0.2, 0.05])
text_box_K = TextBox(axbox_K, "K(x,y) ")
text_box_K.on_submit(submit_fn_K)
text_box_K.set_val(K_str)

def submit_fn_f(value):
    update_function("f_str",str(value))
    plt.draw()
axbox_f = fig.add_axes([0.79, 0.74, 0.2, 0.05])
text_box_f = TextBox(axbox_f, "f(x) ")
text_box_f.on_submit(submit_fn_f)
text_box_f.set_val(f_str)

def submit_fn_ans(value):
    update_function("ans",str(value))
    plt.draw()
axbox_ans = fig.add_axes([0.79, 0.69, 0.2, 0.05])
text_box_ans = TextBox(axbox_ans, "Real ")
text_box_ans.on_submit(submit_fn_ans)
text_box_ans.set_val(ans_str)

plt.show()

