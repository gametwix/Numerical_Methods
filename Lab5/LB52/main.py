import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, TextBox

import PDE
from hyperbolicexplicit import Hyperbolic_explicit_solver as Hyp_exp_solv
from hyperbolicimplicit import Hyperbolic_implicit_solver as Hyp_imp_solv
from nmutils import (
    get_func_x,
    get_func_t,
    get_func_xt,
    count_mse,
    count_mae,
    max_abs_err,
)


class Inteface_params:
    def __init__(self):
        self.start_x = 0
        self.finish_x = np.pi
        self.max_t = 10
        self.resolution_x = 10
        self.resolution_t = 1000

        # U(0, t) ans U(l, t)
        self.Utsx_str = "exp(-t)"
        self.Utfx_str = "-exp(-t)"

        # U(x, 0)
        self.Uxst_str = "sin(x)"
        self.Uxst_str_dt = "-sin(x)"
        self.Uxst_str_dx = "cos(x)"
        self.Uxst_str_dxdx = "-sin(x)"

        self.fi_str = "-cos(x)*exp(-t)"

        self.Uans_str = "exp(-t)*sin(x)"

        self.a = 1
        self.b = 1
        self.c = -1
        self.e = 3
        self.alf0 = 1
        self.bet0 = 0
        self.alfn = 1
        self.betn = 0
        self.lr_method = 1
        self.k1_method = 1
        self.cur_time = 0


# Init interface
p = Inteface_params()

eq = PDE.PartialDifferentialEquation(
    Ut_start_x=get_func_t(p.Utsx_str),
    Ut_finish_x=get_func_t(p.Utfx_str),
    Ux_start_t=get_func_x(p.Uxst_str),
    Ux_start_t_dt=get_func_x(p.Uxst_str_dt),
    Ux_start_t_dx=get_func_x(p.Uxst_str_dx),
    Ux_start_t_dxdx=get_func_x(p.Uxst_str_dxdx),
    fi=get_func_xt(p.fi_str),
    U_answer=get_func_xt(p.Uans_str),
    alf0=p.alf0,
    bet0=p.bet0,
    alfn=p.alfn,
    betn=p.betn,
    a=p.a,
    b=p.b,
    c=p.c,
    e=p.e,
)

explicit = Hyp_exp_solv(
    eq,
    p.start_x,
    p.finish_x,
    p.resolution_x,
    p.max_t,
    p.resolution_t,
    p.lr_method,
    p.k1_method,
)


implicit = Hyp_imp_solv(
    eq,
    p.start_x,
    p.finish_x,
    p.resolution_x,
    p.max_t,
    p.resolution_t,
    p.lr_method,
    p.k1_method,
)

err_t = np.array([i * explicit.tau for i in range(p.resolution_t)])
all_y = [[eq.U_answer(elem, i) for elem in explicit.x] for i in err_t]
p.ind = int(p.cur_time / explicit.tau)
y_ans = all_y[p.ind]


fig = plt.figure("Data")
fig_graph = plt.figure("Graph")
fig_error = plt.figure("Max Absolute Error per Time")
ax_graph = fig_graph.add_subplot()
ax_graph.grid()

ax_error = fig_error.add_subplot()
ax_error.grid()

(real,) = ax_graph.plot(explicit.x, y_ans, label="Аналитическое решение")
(my_plot_exp,) = ax_graph.plot(explicit.x, explicit.answer[p.ind], label="Явная схема")
(my_plot_imp,) = ax_graph.plot(
    implicit.x, implicit.answer[p.ind], label="Неявная схема"
)

exp_max_er = [max_abs_err(all_y[i], explicit.answer[i]) for i in range(len(err_t))]
imp_max_er = [max_abs_err(all_y[i], implicit.answer[i]) for i in range(len(err_t))]
(err_plot_exp,) = ax_error.plot(err_t, exp_max_er, label="Явная схема")
(err_plot_imp,) = ax_error.plot(err_t, imp_max_er, label="Неявная схема")


fig.text(0.75, 0.95, "Errors MSE:")
fig.text(0.75, 0.90, "Explicit method:")
er_exp_mse = fig.text(0.75, 0.87, str(count_mse(explicit.answer[p.ind], y_ans)))
fig.text(0.75, 0.84, "Implicit method:")
er_imp_mse = fig.text(0.75, 0.81, str(count_mse(implicit.answer[p.ind], y_ans)))


fig.text(0.75, 0.75, "Errors MAE:")
fig.text(0.75, 0.70, "Explicit method:")
er_exp_mae = fig.text(0.75, 0.67, str(count_mae(explicit.answer[p.ind], y_ans)))
fig.text(0.75, 0.64, "Implicit method:")
er_imp_mae = fig.text(0.75, 0.61, str(count_mae(implicit.answer[p.ind], y_ans)))

fig.text(0.65, 0.55, "max_t:")
max_t_print = fig.text(0.73, 0.55, str(p.max_t))

fig.text(0.65, 0.52, "σ:")
sig_print = fig.text(0.68, 0.52, str((eq.a * explicit.tau**2) / (explicit.h**2)))


def update_function():
    global eq, p
    global er_exp, er_imp, er_exp_imp, fig_graph, method
    global err_plot_exp, err_plot_imp, err_plot_expimp

    explicit.equation = eq
    explicit.start_x = p.start_x
    explicit.finish_x = p.finish_x
    explicit.resolution_x = p.resolution_x
    explicit.finish_t = p.max_t
    explicit.resolution_t = p.resolution_t
    explicit.lr_method = p.lr_method

    implicit.equation = eq
    implicit.start_x = p.start_x
    implicit.finish_x = p.finish_x
    implicit.resolution_x = p.resolution_x
    implicit.finish_t = p.max_t
    implicit.resolution_t = p.resolution_t
    implicit.lr_method = p.lr_method

    err_t = np.array([i * explicit.tau for i in range(p.resolution_t)])
    all_y = [[eq.U_answer(elem, i) for elem in explicit.x] for i in err_t]
    p.ind = int(p.cur_time / explicit.tau)
    y_ans = all_y[p.ind]
    # print(np.max(ans_exp_imp))
    ans_exp = explicit.answer[p.ind]
    ans_imp = implicit.answer[p.ind]
    real.set_xdata(implicit.x)
    real.set_ydata(y_ans)
    my_plot_exp.set_xdata(explicit.x)
    my_plot_exp.set_ydata(ans_exp)
    my_plot_imp.set_xdata(implicit.x)
    my_plot_imp.set_ydata(ans_imp)

    exp_max_er = [max_abs_err(all_y[i], explicit.answer[i]) for i in range(len(err_t))]
    imp_max_er = [max_abs_err(all_y[i], implicit.answer[i]) for i in range(len(err_t))]
    err_plot_exp.set_xdata(err_t)
    err_plot_exp.set_ydata(exp_max_er)
    err_plot_imp.set_xdata(err_t)
    err_plot_imp.set_ydata(imp_max_er)

    er_exp_mse.set_text(str(count_mse(ans_exp, y_ans)))
    er_imp_mse.set_text(str(count_mse(ans_imp, y_ans)))
    er_exp_mae.set_text(str(count_mae(ans_exp, y_ans)))
    er_imp_mae.set_text(str(count_mae(ans_imp, y_ans)))
    max_t_print.set_text(str(p.max_t))
    # print(explicit.h**2)
    # print(explicit.tau**2)
    sig_print.set_text(str((eq.a * explicit.tau**2) / (explicit.h**2)))

    ax_graph.relim()
    ax_error.relim()
    # ax_graph.set_ylim(0, 1)
    ax_graph.legend()
    ax_error.legend()
    ax_graph.autoscale_view()
    ax_error.autoscale_view()
    fig_graph.canvas.draw()
    fig_error.canvas.draw()


def submit_fn_a(value):
    global eq
    eq.a = float(value)
    update_function()
    plt.draw()


axbox_a = fig.add_axes([0.08, 0.93, 0.05, 0.05])
text_box_a = TextBox(axbox_a, "a ")
text_box_a.on_submit(submit_fn_a)
text_box_a.set_val(eq.a)


def submit_fn_b(value):
    global eq
    eq.b = float(value)
    update_function()
    plt.draw()


axbox_b = fig.add_axes([0.16, 0.93, 0.05, 0.05])
text_box_b = TextBox(axbox_b, "b ")
text_box_b.on_submit(submit_fn_b)
text_box_b.set_val(eq.b)


def submit_fn_c(value):
    global eq
    eq.c = float(value)
    update_function()
    plt.draw()


axbox_c = fig.add_axes([0.24, 0.93, 0.05, 0.05])
text_box_c = TextBox(axbox_c, "c ")
text_box_c.on_submit(submit_fn_c)
text_box_c.set_val(eq.c)


def submit_fn_e(value):
    global eq
    eq.e = float(value)
    update_function()
    plt.draw()


axbox_e = fig.add_axes([0.32, 0.93, 0.05, 0.05])
text_box_e = TextBox(axbox_e, "e ")
text_box_e.on_submit(submit_fn_e)
text_box_e.set_val(eq.e)


def submit_fn_max_x(value):
    global p
    p.finish_x = eval(value)
    update_function()
    plt.draw()


axbox_max_x = fig.add_axes([0.22, 0.51, 0.08, 0.05])
text_box_max_x = TextBox(axbox_max_x, "max_x")
text_box_max_x.on_submit(submit_fn_max_x)
text_box_max_x.set_val("np.pi")


def submit_fn_max_t(value):
    global p
    p.h = float(value)
    p.max_t = p.h * p.resolution_t
    update_function()
    plt.draw()


axbox_max_t = fig.add_axes([0.08, 0.51, 0.07, 0.05])
text_box_max_t = TextBox(axbox_max_t, "tau")
text_box_max_t.on_submit(submit_fn_max_t)
text_box_max_t.set_val(0.01)


def submit_fn_res_x(value):
    global p
    p.resolution_x = int(value)
    update_function()
    plt.draw()


axbox_res_x = fig.add_axes([0.22, 0.44, 0.08, 0.05])
text_box_res_x = TextBox(axbox_res_x, "res_x")
text_box_res_x.on_submit(submit_fn_res_x)
text_box_res_x.set_val(p.resolution_x)


def submit_fn_res_t(value):
    global p
    p.resolution_t = int(value)
    p.max_t = p.h * p.resolution_t
    update_function()
    plt.draw()


axbox_res_t = fig.add_axes([0.08, 0.44, 0.07, 0.05])
text_box_res_t = TextBox(axbox_res_t, "res_t")
text_box_res_t.on_submit(submit_fn_res_t)
text_box_res_t.set_val(p.resolution_t)


def submit_fn_fi(value):
    global eq
    eq._f_fi = get_func_xt(value)
    update_function()
    plt.draw()


axbox_fi = fig.add_axes([0.08, 0.86, 0.21, 0.05])
text_box_fi = TextBox(axbox_fi, "f(x,t) ")
text_box_fi.on_submit(submit_fn_fi)
text_box_fi.set_val(p.fi_str)


def submit_fn_uxst(value):
    global p, eq
    p.Uxst_str = value
    eq.Ux_start_t = get_func_x(value)
    update_function()
    plt.draw()


axbox_uxst = fig.add_axes([0.08, 0.79, 0.21, 0.05])
text_box_uxst = TextBox(axbox_uxst, "U(x,0)")
text_box_uxst.on_submit(submit_fn_uxst)
text_box_uxst.set_val(p.Uxst_str)


def submit_fn_uxstdt(value):
    global p, eq
    p.Uxst_str_dt = value
    eq.Ux_start_t_dt = get_func_x(value)
    update_function()
    plt.draw()


axbox_uxstdt = fig.add_axes([0.40, 0.79, 0.21, 0.05])
text_box_uxstdt = TextBox(axbox_uxstdt, "U(x,0)ₜ")
text_box_uxstdt.on_submit(submit_fn_uxstdt)
text_box_uxstdt.set_val(p.Uxst_str_dt)


def submit_fn_uxstdx(value):
    global p, eq
    p.Uxst_str_dx = value
    eq.Ux_start_t_dx = get_func_x(value)
    update_function()
    plt.draw()


axbox_uxstdx = fig.add_axes([0.40, 0.72, 0.21, 0.05])
text_box_uxstdx = TextBox(axbox_uxstdx, "U(x,0)ₓ")
text_box_uxstdx.on_submit(submit_fn_uxstdx)
text_box_uxstdx.set_val(p.Uxst_str_dx)


def submit_fn_uxstdxdx(value):
    global p, eq
    p.Uxst_str_dxdx = value
    eq.Ux_start_t_dxdx = get_func_x(value)
    update_function()
    plt.draw()


axbox_uxstdxdx = fig.add_axes([0.40, 0.64, 0.21, 0.05])
text_box_uxstdxdx = TextBox(axbox_uxstdxdx, "U(x,0)ₓₓ")
text_box_uxstdxdx.on_submit(submit_fn_uxstdxdx)
text_box_uxstdxdx.set_val(p.Uxst_str_dxdx)


def submit_fn_utsx(value):
    global p, eq
    p.Utsx_str = value
    eq._f_Ut_start_x = get_func_t(value)
    update_function()
    plt.draw()


axbox_utsx = fig.add_axes([0.08, 0.72, 0.21, 0.05])
text_box_utsx = TextBox(axbox_utsx, "φ(0,t)")
text_box_utsx.on_submit(submit_fn_utsx)
text_box_utsx.set_val(p.Utsx_str)


def submit_fn_utfx(value):
    global p, eq
    p.Utfx_str = value
    eq._f_Ut_finish_x = get_func_t(value)
    update_function()
    plt.draw()


axbox_utfx = fig.add_axes([0.08, 0.65, 0.21, 0.05])
text_box_utfx = TextBox(axbox_utfx, "φ(l,t)")
text_box_utfx.on_submit(submit_fn_utfx)
text_box_utfx.set_val(p.Utfx_str)


def submit_fn_ans(value):
    global p, eq
    p.Uans_str = value
    eq._f_U_answer = get_func_xt(value)
    update_function()
    plt.draw()


axbox_ans = fig.add_axes([0.08, 0.58, 0.21, 0.05])
text_box_ans = TextBox(axbox_ans, "U_ans")
text_box_ans.on_submit(submit_fn_ans)
text_box_ans.set_val(p.Uans_str)


def submit_fn_slider(value):
    global p
    p.cur_time = p.max_t * float(value)
    update_function()
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
    global eq, p
    p.alf0 = float(value)
    eq.alf0 = float(value)
    update_function()
    plt.draw()


axbox_alph = fig.add_axes([0.08, 0.37, 0.05, 0.05])
text_box_alph = TextBox(axbox_alph, "α ")
text_box_alph.on_submit(submit_fn_alph)
text_box_alph.set_val(eq.alf0)


def submit_fn_bet(value):
    global eq, p
    p.bet0 = float(value)
    eq.bet0 = float(value)
    update_function()
    plt.draw()


axbox_bet = fig.add_axes([0.16, 0.37, 0.05, 0.05])
text_box_bet = TextBox(axbox_bet, "β ")
text_box_bet.on_submit(submit_fn_bet)
text_box_bet.set_val(eq.bet0)


def submit_fn_gam(value):
    global eq, p
    p.alfn = float(value)
    eq.alfn = float(value)
    update_function()
    plt.draw()


axbox_gam = fig.add_axes([0.08, 0.3, 0.05, 0.05])
text_box_gam = TextBox(axbox_gam, "γ ")
text_box_gam.on_submit(submit_fn_gam)
text_box_gam.set_val(eq.alfn)


def submit_fn_delt(value):
    global eq, p
    p.betn = float(value)
    eq.betn = float(value)
    update_function()
    plt.draw()


axbox_delt = fig.add_axes([0.16, 0.3, 0.05, 0.05])
text_box_delt = TextBox(axbox_delt, "δ ")
text_box_delt.on_submit(submit_fn_delt)
text_box_delt.set_val(eq.betn)


def submit_fn_met_lr(value):
    global p
    p.lr_method = int(value)
    update_function()
    plt.draw()


axbox_met_lr = fig.add_axes([0.08, 0.23, 0.05, 0.05])
text_box_met_lr = TextBox(axbox_met_lr, "metlr ")
text_box_met_lr.on_submit(submit_fn_met_lr)
text_box_met_lr.set_val(p.lr_method)


def submit_fn_met_k1(value):
    global p
    p.k1_method = int(value)
    update_function()
    plt.draw()


axbox_met_k1 = fig.add_axes([0.22, 0.23, 0.05, 0.05])
text_box_met_k1 = TextBox(axbox_met_k1, "metk1 ")
text_box_met_k1.on_submit(submit_fn_met_k1)
text_box_met_k1.set_val(p.k1_method)
plt.show()
