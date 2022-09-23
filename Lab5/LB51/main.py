from cProfile import label
from math import sin, cos, exp

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox 
import numpy as np

import PDE


def get_func_x(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")
    def func(x, **kwargs):
        kwargs['x'] = x
        return eval(code, {'sin': sin, 'cos': cos, 'exp': exp}, kwargs)
    return func


def get_func_t(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")
    def func(t, **kwargs):
        kwargs['t'] = t
        return eval(code, {'sin': sin, 'cos': cos, 'exp': exp}, kwargs)
    return func


def get_func_xt(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")
    def func(x, t, **kwargs):
        kwargs['x'] = x
        kwargs['t'] = t
        
        return eval(code, {'sin': sin, 'cos': cos, 'exp': exp}, kwargs)
    return func


def count_MSE(y,y_ans):
    ans = 0
    for i in range(len(y)):
        ans += (y[i] - y_ans[i])**2
        
    return np.sqrt(ans)

start_x = 0
finish_x = np.pi
max_t = 1
resolution_x = 30
resolution_t = 10

Uxst_str = 'sin(x)'
Uans_str = 'exp(-a*t)*sin(x)'
Utsx_str = '0'
Utfx_str = '0'
fi_str = '0'




eq = PDE.PartialDifferentialEquation(Ux_start_t=get_func_x(Uxst_str),
                                     U_answer=get_func_xt(Uans_str),
                                     Ut_start_x=get_func_t(Utsx_str),
                                     Ut_finish_x=get_func_t(Utfx_str),
                                     fi=get_func_xt(fi_str)
                                     )


h = (finish_x - start_x) / (resolution_x - 1)
tau = max_t / (resolution_t - 1)

x = np.array([start_x + h*i for i in range(resolution_x)])
y_ans = [eq.U_answer(elem, 0) for elem in x]
cur_time = 0


expimplicit = PDE.Explicit_Implicit_method(eq, start_x, finish_x,
                                        resolution_x, max_t, resolution_t)
explicit = PDE.Explicit_method(eq, start_x, finish_x,
                                        resolution_x, max_t, resolution_t)
implicit = PDE.Implicit_method(eq, start_x, finish_x,
                                        resolution_x, max_t, resolution_t)

fig = plt.figure()
fig.subplots_adjust(right=0.7, bottom=0.2, top=0.95, left=0.05)
ax = fig.add_subplot()
ax.grid()
real, = ax.plot(x, y_ans, label='Аналитическое решение')
my_plot_exp_imp, = ax.plot(x, expimplicit[0], label='Кранк-Николсона')
my_plot_exp, = ax.plot(x, explicit[0], label='Явная схема')
my_plot_imp, = ax.plot(x, implicit[0], label='Неявная схема')



fig.text(0.84,0.4, 'Errors:')
fig.text(0.72,0.35, 'Explicit method:')
er_exp = fig.text(0.72,0.32, str(count_MSE(explicit[0], y_ans)))
fig.text(0.72,0.29, 'Implicit method:')
er_imp = fig.text(0.72,0.26, str(count_MSE(implicit[0], y_ans)))

fig.text(0.72,0.23, 'Crank-Nicolson method:')
er_exp_imp = fig.text(0.72,0.2, str(count_MSE(expimplicit[0], y_ans)))

def update_function(param_name, param_val):
    global eq, x, cur_time, finish_x, max_t, resolution_x, resolution_t, tau, h,er_exp, er_imp, er_exp_imp
    if param_name == 'a':
        eq.a = param_val
    elif param_name == 'b':
        eq.b = param_val
    elif param_name == 'c':
        eq.c = param_val
    elif param_name == 'slider':
        cur_time = max_t*param_val
    elif param_name == 'fi':
        fi_str = param_val
        eq._f_fi = get_func_xt(fi_str)
    elif param_name == 'uxst':
        Uxst_str = param_val
        eq._f_Ux_start_t = get_func_x(Uxst_str)
    elif param_name == 'utsx':
        Utsx_str = param_val
        eq._f_Ut_start_x = get_func_t(Utsx_str)
    elif param_name == 'utfx':
        Utfx_str = param_val
        eq._f_Ut_finish_x = get_func_t(Utfx_str)
    elif param_name == 'uans':
        Uans_str = param_val
        eq._f_U_answer = get_func_xt(Uans_str)
    elif param_name == 'max_x':
        finish_x = eval(param_val)
    elif param_name == 'max_t':
        max_t = param_val

    elif param_name == 'res_x':
        resolution_x = param_val
    elif param_name == 'res_t':
        resolution_t = param_val
        
    h = (finish_x - start_x) / (resolution_x - 1)
    tau = max_t / (resolution_t - 1)
    x = np.array([start_x + h*i for i in range(resolution_x)])
    y_ans = [eq.U_answer(elem, int(cur_time / tau) * tau) for elem in x]
    
    expimplicit = PDE.Explicit_Implicit_method(eq, start_x, finish_x,
                                            resolution_x, max_t, resolution_t)
    explicit = PDE.Explicit_method(eq, start_x, finish_x,
                                        resolution_x, max_t, resolution_t)
    implicit = PDE.Implicit_method(eq, start_x, finish_x,
                                        resolution_x, max_t, resolution_t)
    
    ans_exp_imp = expimplicit[int(cur_time / tau)]
    ans_exp = explicit[int(cur_time / tau)]
    ans_imp = implicit[int(cur_time / tau)]
    real.set_xdata(x)
    real.set_ydata(y_ans)
    my_plot_exp_imp.set_xdata(x)
    my_plot_exp_imp.set_ydata(ans_exp_imp)
    my_plot_exp.set_xdata(x)
    my_plot_exp.set_ydata(ans_exp)
    my_plot_imp.set_xdata(x)
    my_plot_imp.set_ydata(ans_imp)
    
    er_exp.set_text(str(count_MSE(ans_exp, y_ans)))
    er_imp.set_text(str(count_MSE(ans_imp, y_ans)))
    er_exp_imp.set_text(str(count_MSE(ans_exp_imp, y_ans)))
    
    ax.relim()
    ax.set_ylim(0,1)
    ax.legend()
    ax.autoscale_view()


def submit_fn_a(value):
    update_function('a', float(value))
    plt.draw()


axbox_a = fig.add_axes([0.77, 0.93, 0.05, 0.05])
text_box_a = TextBox(axbox_a, "a ")
text_box_a.on_submit(submit_fn_a)
text_box_a.set_val(eq.a)

def submit_fn_b(value):
    update_function('b', float(value))
    plt.draw()


axbox_b = fig.add_axes([0.85, 0.93, 0.05, 0.05])
text_box_b = TextBox(axbox_b, "b ")
text_box_b.on_submit(submit_fn_b)
text_box_b.set_val(eq.b)

def submit_fn_c(value):
    update_function('c', float(value))
    plt.draw()
    

axbox_c = fig.add_axes([0.93, 0.93, 0.05, 0.05])
text_box_c = TextBox(axbox_c, "c ")
text_box_c.on_submit(submit_fn_c)
text_box_c.set_val(eq.c)


def submit_fn_max_x(value):
    update_function('max_x', value)
    plt.draw()


axbox_max_x = fig.add_axes([0.77, 0.51, 0.08, 0.05])
text_box_max_x = TextBox(axbox_max_x, "max_x")
text_box_max_x.on_submit(submit_fn_max_x)
text_box_max_x.set_val('np.pi')


def submit_fn_max_t(value):
    update_function('max_t', float(value))
    plt.draw()


axbox_max_t = fig.add_axes([0.93, 0.51, 0.05, 0.05])
text_box_max_t = TextBox(axbox_max_t, "max_t")
text_box_max_t.on_submit(submit_fn_max_t)
text_box_max_t.set_val(1)

def submit_fn_res_x(value):
    update_function('res_x', int(value))
    plt.draw()


axbox_res_x = fig.add_axes([0.77, 0.44, 0.08, 0.05])
text_box_res_x = TextBox(axbox_res_x, "res_x")
text_box_res_x.on_submit(submit_fn_res_x)
text_box_res_x.set_val(resolution_x)


def submit_fn_res_t(value):
    update_function('res_t', int(value))
    plt.draw()


axbox_res_t = fig.add_axes([0.91, 0.44, 0.07, 0.05])
text_box_res_t = TextBox(axbox_res_t, "res_t")
text_box_res_t.on_submit(submit_fn_res_t)
text_box_res_t.set_val(resolution_t)


def submit_fn_fi(value):
    update_function('fi', value)
    plt.draw()


axbox_fi = fig.add_axes([0.77, 0.86, 0.21, 0.05])
text_box_fi = TextBox(axbox_fi, "f(x,t) ")
text_box_fi.on_submit(submit_fn_fi)
text_box_fi.set_val(fi_str)

def submit_fn_uxst(value):
    update_function('uxst', value)
    plt.draw()


axbox_uxst = fig.add_axes([0.77, 0.79, 0.21, 0.05])
text_box_uxst = TextBox(axbox_uxst, "U(x,0)")
text_box_uxst.on_submit(submit_fn_uxst)
text_box_uxst.set_val(Uxst_str)

def submit_fn_utsx(value):
    update_function('utsx', value)
    plt.draw()


axbox_utsx = fig.add_axes([0.77, 0.72, 0.21, 0.05])
text_box_utsx = TextBox(axbox_utsx, "U(0,t)")
text_box_utsx.on_submit(submit_fn_utsx)
text_box_utsx.set_val(Utsx_str)

def submit_fn_utfx(value):
    update_function('utfx', value)
    plt.draw()


axbox_utfx = fig.add_axes([0.77, 0.65, 0.21, 0.05])
text_box_utfx = TextBox(axbox_utfx, "U(l,t)")
text_box_utfx.on_submit(submit_fn_utfx)
text_box_utfx.set_val(Utfx_str)

def submit_fn_ans(value):
    update_function('uans', value)
    plt.draw()


axbox_ans = fig.add_axes([0.77, 0.58, 0.21, 0.05])
text_box_ans = TextBox(axbox_ans, "U_ans")
text_box_ans.on_submit(submit_fn_ans)
text_box_ans.set_val(Uans_str)


def submit_fn_slider(value):
    update_function('slider', float(value))
    plt.draw()

axbox_time = fig.add_axes([0.1, 0.05, 0.8, 0.05])
time_slider = Slider(
    ax=axbox_time,
    label='time',
    valmin=0,
    valmax=1,
    valinit=0,
)



time_slider.on_changed(submit_fn_slider)


plt.show()

