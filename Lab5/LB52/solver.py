import numpy as np


class Solver:
    def __init__(
        self,
        equation,
        start_x,
        finish_x,
        resolution_x,
        finish_t,
        resolution_t,
        lr_method,
        k1_method,
    ):
        self._equation = equation
        self._start_x = start_x
        self._finish_x = finish_x
        self._resolution_x = resolution_x
        self._finish_t = finish_t
        self._resolution_t = resolution_t
        self._lr_method = lr_method
        self._k1_method = k1_method
        self._h = self._calc_h()
        self._tau = self._calc_tau()
        self._x = self._calc_x()

    @property
    def equation(self):
        return self._equation

    @equation.setter
    def equation(self, value):
        self._equation = value
        if hasattr(self, "_answer"):
            delattr(self, "_answer")

    @property
    def start_x(self):
        return self._start_x

    @start_x.setter
    def start_x(self, value):
        if self._start_x != value:
            self._start_x = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")
        
    @property
    def finish_x(self):
        return self._finish_x

    @finish_x.setter
    def finish_x(self, value):
        if self._finish_x != value:
            self._finish_x = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")
    
    @property
    def resolution_x(self):
        return self._resolution_x

    @resolution_x.setter
    def resolution_x(self, value):
        if self._resolution_x != value:
            self._resolution_x = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")

    @property
    def finish_t(self):
        return self._finish_t

    @finish_t.setter
    def finish_t(self, value):
        if self._finish_t != value:
            self._finish_t = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")

    @property
    def resolution_t(self):
        return self._resolution_t

    @resolution_t.setter
    def resolution_t(self, value):
        if self._resolution_t != value:
            self._resolution_t = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")

    @property
    def lr_method(self):
        return self._lr_method

    @lr_method.setter
    def lr_method(self, value):
        if self._lr_method != value:
            self._lr_method = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")

    @property
    def k1_method(self):
        return self._k1_method

    @k1_method.setter
    def k1_method(self, value):
        if self._k1_method != value:
            self._k1_method = value
            if hasattr(self, "_answer"):
                delattr(self, "_answer")

    @property
    def h(self):
        return self._h

    @property
    def tau(self):
        return self._tau

    @property
    def x(self):
        return self._x

    @property
    def answer(self):
        self._h = self._calc_h()
        self._tau = self._calc_tau()
        self._x = self._calc_x()
        if not hasattr(self, "_answer"):
            self._calc_answer()
        return self._answer

    def _calc_h(self):
        return (self._finish_x - self._start_x) / (self._resolution_x - 1)

    def _calc_tau(self):
        return self.finish_t / (self.resolution_t - 1)

    def _calc_x(self):
        return np.array([self._start_x + self._h * i for i in range(self._resolution_x)])

    def _calc_answer(self):
        pass
