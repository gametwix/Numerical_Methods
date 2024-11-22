import numpy as np
from math import sin, cos, exp


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


def count_mse(y, y_ans):
    return np.sqrt(np.sum((y - y_ans) ** 2))


def count_mae(y, y_ans):
    return np.sum(np.abs(y - y_ans))


def max_abs_err(y, y_ans):
    return np.max(np.abs(y - y_ans))
