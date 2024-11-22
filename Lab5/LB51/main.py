from math import sin, cos, exp
from os import EX_IOERR

import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox
import numpy as np

import PDE


def get_func_x(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")

    def func(x, **kwargs):
        kwargs["x"] = x
        return eval(code, {"sin": sin, "cos": cos, "exp": exp}, kwargs)

    return func


def get_func_t(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")

    def func(t, **kwargs):
        kwargs["t"] = t
        return eval(code, {"sin": sin, "cos": cos, "exp": exp}, kwargs)

    return func


def get_func_xt(str_math_func: str):
    code = compile(str_math_func, "<string>", "eval")

    def func(x, t, **kwargs):
        kwargs["x"] = x
        kwargs["t"] = t

        return eval(code, {"sin": sin, "cos": cos, "exp": exp}, kwargs)

    return func


def count_MSE(y, y_ans):
    ans = 0
    for i in range(len(y)):
        ans += (y[i] - y_ans[i]) ** 2

    return np.sqrt(ans)


def count_MAE(y, y_ans):
    ans = 0
    for i in range(len(y)):
        ans += np.abs(y[i] - y_ans[i]) 

    return np.abs(ans)

def MaxAbsErr(y, y_ans):
    m = 0
    for i in range(len(y)):
        if np.abs(y[i] - y_ans[i]) > m:
            m = np.abs(y[i] - y_ans[i])
    return m


start_x = 0
finish_x = np.pi
max_t = 1
resolution_x = 30
resolution_t = 200

Uxst_str = "sin(x)"
Uans_str = "exp(-a*t)*sin(x)"
Utsx_str = "exp(-a*t)"
Utfx_str = "-exp(-a*t)"
# Utsx_str = "0"
# Utfx_str = "0"
fi_str = "0"


eq = PDE.PartialDifferentialEquation(
    Ux_start_t=get_func_x(Uxst_str),
    U_answer=get_func_xt(Uans_str),
    Ut_start_x=get_func_t(Utsx_str),
    Ut_finish_x=get_func_t(Utfx_str),
    fi=get_func_xt(fi_str),
    alf0 = 1,
    bet0 = 0,
    alfn= 1,
    betn=0
)


h = (finish_x - start_x) / (resolution_x - 1)
tau = max_t / (resolution_t - 1)

x = np.array([start_x + h * i for i in range(resolution_x)])
y_ans = [eq.U_answer(elem, 0) for elem in x]

cur_time = 0

method = 3


expimplicit = PDE.Explicit_Implicit_method(
    eq, start_x, finish_x, resolution_x, max_t, resolution_t,method=method
)
explicit = PDE.Explicit_method(eq,
                               start_x,
                               finish_x,
                               resolution_x,
                               max_t,
                               resolution_t,
                               method=method)
implicit = PDE.Implicit_method(eq, start_x,
                               finish_x,
                               resolution_x,
                               max_t,
                               resolution_t,
                               method=method)

fig = plt.figure('Data')
fig_graph = plt.figure('Graph')
fig_error = plt.figure('Max Absolute Error per Time')
ax_graph = fig_graph.add_subplot()
ax_graph.grid()

ax_error = fig_error.add_subplot()
ax_error.grid()

(real,) = ax_graph.plot(x, y_ans, label="Аналитическое решение")
(my_plot_exp_imp,) = ax_graph.plot(x, expimplicit[0], label="Кранк-Николсона")
(my_plot_exp,) = ax_graph.plot(x, explicit[0], label="Явная схема")
(my_plot_imp,) = ax_graph.plot(x, implicit[0], label="Неявная схема")


err_t = np.array([i*tau for i in range(resolution_t)])
all_y = [[eq.U_answer(elem, i) for elem in x] for i in err_t]

exp_max_er = [MaxAbsErr(all_y[i], explicit[i]) for i in range(len(err_t))]
imp_max_er = [MaxAbsErr(all_y[i], implicit[i]) for i in range(len(err_t))]
expimp_max_er = [MaxAbsErr(all_y[i], expimplicit[i]) for i in range(len(err_t))]
(err_plot_exp,) = ax_error.plot(err_t, exp_max_er, label="Явная схема")
(err_plot_imp,) = ax_error.plot(err_t, imp_max_er, label="Неявная схема")
(err_plot_exp_imp,) = ax_error.plot(err_t, expimp_max_er, label="Кранк-Николсона")






fig.text(0.4, 0.95, "Errors MSE:")
fig.text(0.4, 0.90, "Explicit method:")
er_exp_mse = fig.text(0.4, 0.87, str(count_MSE(explicit[0], y_ans)))
fig.text(0.4, 0.84, "Implicit method:")
er_imp_mse = fig.text(0.4, 0.81, str(count_MSE(implicit[0], y_ans)))

fig.text(0.4, 0.78, "Crank-Nicolson method:")
er_exp_imp_mse = fig.text(0.4, 0.75, str(count_MSE(expimplicit[0], y_ans)))

fig.text(0.4, 0.55, "Errors MAE:")
fig.text(0.4, 0.50, "Explicit method:")
er_exp_mae = fig.text(0.4, 0.47, str(count_MAE(explicit[0], y_ans)))
fig.text(0.4, 0.44, "Implicit method:")
er_imp_mae = fig.text(0.4, 0.41, str(count_MAE(implicit[0], y_ans)))

fig.text(0.4, 0.38, "Crank-Nicolson method:")
er_exp_imp_mae = fig.text(0.4, 0.35, str(count_MAE(expimplicit[0], y_ans)))


def update_function(param_name, param_val):
    global eq, x, cur_time, finish_x, max_t, resolution_x, resolution_t, tau, h
    global er_exp, er_imp, er_exp_imp, fig_graph, method
    global err_plot_exp, err_plot_imp, err_plot_expimp
    if param_name == "a":
        eq.a = param_val
    elif param_name == "b":
        eq.b = param_val
    elif param_name == "c":
        eq.c = param_val
    elif param_name == "slider":
        cur_time = max_t * param_val
    elif param_name == "fi":
        fi_str = param_val
        eq._f_fi = get_func_xt(fi_str)
    elif param_name == "uxst":
        Uxst_str = param_val
        eq._f_Ux_start_t = get_func_x(Uxst_str)
    elif param_name == "utsx":
        Utsx_str = param_val
        eq._f_Ut_start_x = get_func_t(Utsx_str)
    elif param_name == "utfx":
        Utfx_str = param_val
        eq._f_Ut_finish_x = get_func_t(Utfx_str)
    elif param_name == "uans":
        Uans_str = param_val
        eq._f_U_answer = get_func_xt(Uans_str)
    elif param_name == "max_x":
        finish_x = eval(param_val)
    elif param_name == "max_t":
        max_t = param_val
    elif param_name == "res_x":
        resolution_x = param_val
    elif param_name == "res_t":
        resolution_t = param_val
    elif param_name == "alph":
        eq.alf0 = param_val
    elif param_name == "bet":
        eq.bet0 = param_val
    elif param_name == "gam":
        eq.alfn = param_val
    elif param_name == "delt":
        eq.betn = param_val
    elif param_name == "met":
        method = param_val

    h = (finish_x - start_x) / (resolution_x - 1)
    tau = max_t / (resolution_t - 1)
    x = np.array([start_x + h * i for i in range(resolution_x)])
    y_ans = [eq.U_answer(elem, int(cur_time / tau) * tau) for elem in x]

    expimplicit = PDE.Explicit_Implicit_method(
        eq, start_x, finish_x, resolution_x, max_t, resolution_t, method=method
    )
    explicit = PDE.Explicit_method(
        eq, start_x, finish_x, resolution_x, max_t, resolution_t, method=method
    )
    implicit = PDE.Implicit_method(
        eq, start_x, finish_x, resolution_x, max_t, resolution_t, method=method
    )

    
    ans_exp_imp = expimplicit[int(cur_time / tau)]
    # print(np.max(ans_exp_imp))
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
    
    err_t = np.array([i*tau for i in range(resolution_t)])
    all_y = [[eq.U_answer(elem, i) for elem in x] for i in err_t]

    exp_max_er = [MaxAbsErr(all_y[i], explicit[i]) for i in range(len(err_t))]
    imp_max_er = [MaxAbsErr(all_y[i], implicit[i]) for i in range(len(err_t))]
    expimp_max_er = [MaxAbsErr(all_y[i], expimplicit[i]) for i in range(len(err_t))]
    err_plot_exp.set_xdata(err_t)
    err_plot_exp.set_ydata(exp_max_er)
    err_plot_imp.set_xdata(err_t)
    err_plot_imp.set_ydata(imp_max_er)
    err_plot_exp_imp.set_xdata(err_t)
    err_plot_exp_imp.set_ydata(expimp_max_er)
    
    er_exp_mse.set_text(str(count_MSE(ans_exp, y_ans)))
    er_imp_mse.set_text(str(count_MSE(ans_imp, y_ans)))
    er_exp_imp_mse.set_text(str(count_MSE(ans_exp_imp, y_ans)))
    er_exp_mae.set_text(str(count_MAE(ans_exp, y_ans)))
    er_imp_mae.set_text(str(count_MAE(ans_imp, y_ans)))
    er_exp_imp_mae.set_text(str(count_MAE(ans_exp_imp, y_ans)))

    ax_graph.relim()
    ax_error.relim()
    ax_graph.set_ylim(0, 1)
    ax_graph.legend()
    ax_error.legend()
    ax_graph.autoscale_view()
    ax_error.autoscale_view()
    fig_graph.canvas.draw()
    fig_error.canvas.draw()
    
    
    



def submit_fn_a(value):
    update_function("a", float(value))
    plt.draw()


axbox_a = fig.add_axes([0.08, 0.93, 0.05, 0.05])
text_box_a = TextBox(axbox_a, "a ")
text_box_a.on_submit(submit_fn_a)
text_box_a.set_val(eq.a)


def submit_fn_b(value):
    update_function("b", float(value))
    plt.draw()


axbox_b = fig.add_axes([0.16, 0.93, 0.05, 0.05])
text_box_b = TextBox(axbox_b, "b ")
text_box_b.on_submit(submit_fn_b)
text_box_b.set_val(eq.b)


def submit_fn_c(value):
    update_function("c", float(value))
    plt.draw()


axbox_c = fig.add_axes([0.24, 0.93, 0.05, 0.05])
text_box_c = TextBox(axbox_c, "c ")
text_box_c.on_submit(submit_fn_c)
text_box_c.set_val(eq.c)


def submit_fn_max_x(value):
    update_function("max_x", value)
    plt.draw()


axbox_max_x = fig.add_axes([0.22, 0.51, 0.08, 0.05])
text_box_max_x = TextBox(axbox_max_x, "max_x")
text_box_max_x.on_submit(submit_fn_max_x)
text_box_max_x.set_val("np.pi")


def submit_fn_max_t(value):
    update_function("max_t", float(value))
    plt.draw()


axbox_max_t = fig.add_axes([0.08, 0.51, 0.05, 0.05])
text_box_max_t = TextBox(axbox_max_t, "max_t")
text_box_max_t.on_submit(submit_fn_max_t)
text_box_max_t.set_val(1)


def submit_fn_res_x(value):
    update_function("res_x", int(value))
    plt.draw()


axbox_res_x = fig.add_axes([0.22, 0.44, 0.08, 0.05])
text_box_res_x = TextBox(axbox_res_x, "res_x")
text_box_res_x.on_submit(submit_fn_res_x)
text_box_res_x.set_val(resolution_x)


def submit_fn_res_t(value):
    update_function("res_t", int(value))
    plt.draw()


axbox_res_t = fig.add_axes([0.08, 0.44, 0.07, 0.05])
text_box_res_t = TextBox(axbox_res_t, "res_t")
text_box_res_t.on_submit(submit_fn_res_t)
text_box_res_t.set_val(resolution_t)


def submit_fn_fi(value):
    update_function("fi", value)
    plt.draw()


axbox_fi = fig.add_axes([0.08, 0.86, 0.21, 0.05])
text_box_fi = TextBox(axbox_fi, "f(x,t) ")
text_box_fi.on_submit(submit_fn_fi)
text_box_fi.set_val(fi_str)


def submit_fn_uxst(value):
    update_function("uxst", value)
    plt.draw()


axbox_uxst = fig.add_axes([0.08, 0.79, 0.21, 0.05])
text_box_uxst = TextBox(axbox_uxst, "U(x,0)")
text_box_uxst.on_submit(submit_fn_uxst)
text_box_uxst.set_val(Uxst_str)


def submit_fn_utsx(value):
    update_function("utsx", value)
    plt.draw()


axbox_utsx = fig.add_axes([0.08, 0.72, 0.21, 0.05])
text_box_utsx = TextBox(axbox_utsx, "φ(0,t)")
text_box_utsx.on_submit(submit_fn_utsx)
text_box_utsx.set_val(Utsx_str)


def submit_fn_utfx(value):
    update_function("utfx", value)
    plt.draw()


axbox_utfx = fig.add_axes([0.08, 0.65, 0.21, 0.05])
text_box_utfx = TextBox(axbox_utfx, "φ(l,t)")
text_box_utfx.on_submit(submit_fn_utfx)
text_box_utfx.set_val(Utfx_str)


def submit_fn_ans(value):
    update_function("uans", value)
    plt.draw()


axbox_ans = fig.add_axes([0.08, 0.58, 0.21, 0.05])
text_box_ans = TextBox(axbox_ans, "U_ans")
text_box_ans.on_submit(submit_fn_ans)
text_box_ans.set_val(Uans_str)


def submit_fn_slider(value):
    update_function("slider", float(value))
    plt.draw()


axbox_time = fig.add_axes([0.1, 0.05, 0.8, 0.05])
time_slider = Slider(
    ax=axbox_time,
    label="time",
    valmin=0,
    valmax=1,
    valinit=0,
)


time_slider.on_changed(submit_fn_slider)




def submit_fn_alph(value):
    update_function("alph", float(value))
    plt.draw()


axbox_alph = fig.add_axes([0.08, 0.37, 0.05, 0.05])
text_box_alph = TextBox(axbox_alph, "α ")
text_box_alph.on_submit(submit_fn_alph)
text_box_alph.set_val(eq.alf0)


def submit_fn_bet(value):
    update_function("bet", float(value))
    plt.draw()


axbox_bet = fig.add_axes([0.16, 0.37, 0.05, 0.05])
text_box_bet = TextBox(axbox_bet, "β ")
text_box_bet.on_submit(submit_fn_bet)
text_box_bet.set_val(eq.bet0)


def submit_fn_gam(value):
    update_function("gam", float(value))
    plt.draw()


axbox_gam = fig.add_axes([0.08, 0.3, 0.05, 0.05])
text_box_gam = TextBox(axbox_gam, "γ ")
text_box_gam.on_submit(submit_fn_gam)
text_box_gam.set_val(eq.alfn)


def submit_fn_delt(value):
    update_function("delt", float(value))
    plt.draw()


axbox_delt = fig.add_axes([0.16, 0.3, 0.05, 0.05])
text_box_delt = TextBox(axbox_delt, "δ ")
text_box_delt.on_submit(submit_fn_delt)
text_box_delt.set_val(eq.betn)


def submit_fn_met(value):
    update_function("met", int(value))
    plt.draw()


axbox_met = fig.add_axes([0.08, 0.23, 0.05, 0.05])
text_box_met = TextBox(axbox_met, "met ")
text_box_met.on_submit(submit_fn_met)
text_box_met.set_val(method)


plt.show()
