import numpy as np
import solver


class Hyperbolic_explicit_solver(solver.Solver):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _k1_method_two_point(self):
        if not hasattr(self, "_answer"):
            raise RuntimeError
        for i in range(self._resolution_x):
            self._answer[1][i] = (
                self._answer[0][i]
                + self._equation.Ux_start_t_dt(self._x[i]) * self._tau
            )

    def _k1_method_taylor_series(self):
        if not hasattr(self, "_answer"):
            raise RuntimeError
        for i in range(self._resolution_x):
            self._answer[1][i] = (
                self._answer[0][i]
                + self._equation.Ux_start_t_dt(self._x[i]) * self._tau
                + self._tau**2
                / 2
                * (
                    self._equation.a * self._equation.Ux_start_t_dxdx(self._x[i])
                    + self._equation.b * self._equation.Ux_start_t_dx(self._x[i])
                    + self._equation.c * self._equation.Ux_start_t(self._x[i])
                    + self._equation.fi(self._x[i], 0)
                    - self._equation.e * self._equation.Ux_start_t_dt(self._x[i])
                )
            )

    def _calc_midle_points(self, k):
        if not hasattr(self, "_answer"):
            raise RuntimeError

        denominator = self._tau * self._equation.e * self._h**2 + self._h**2
        fraction_h = self._h**2 / denominator
        fraction_a = self._tau**2 * self._equation.a / denominator
        fraction_b = self._tau**2 * self._equation.b * self._h / denominator
        fraction_c = self._tau**2 * self._equation.c * fraction_h
        fraction_e = self._tau * self._equation.e * fraction_h

        for i in range(1, self._resolution_x - 1):
            self._answer[k + 1][i] = (
                (2 * fraction_h + fraction_c + fraction_e - 2 * fraction_a)
                * self._answer[k][i]
                + (fraction_a + fraction_b / 2) * self._answer[k][i + 1]
                + (fraction_a - fraction_b / 2) * self._answer[k][i - 1]
                - self._answer[k - 1][i] * fraction_h
                + self._tau**2
                * fraction_h
                * self.equation.fi(self._x[i], self._tau * k)
            )

    def _two_points_first_order(self, k):
        denominator_l = self._equation.alf0 - self._equation.bet0 * self._h
        denominator_r = self._equation.betn * self._h + self._equation.alfn

        self._answer[k + 1][0] = (
            self._answer[k + 1][1] * self._equation.alf0
            - self._equation.Ut_start_x(self._tau * (k + 1)) * self._h
        ) / denominator_l

        self._answer[k + 1][-1] = (
            self._equation.Ut_finish_x(self._tau * (k + 1)) * self._h
            + self._answer[k + 1][-2] * self._equation.alfn
        ) / denominator_r

    def _three_points_second_order(self, k):
        denominator_l = 3 * self._equation.alf0 - 2 * self._equation.bet0 * self._h
        denominator_r = 3 * self._equation.alfn + 2 * self._equation.betn * self._h

        self._answer[k + 1][0] = (
            4 * self._answer[k + 1][1] * self._equation.alf0
            - self._answer[k + 1][2] * self._equation.alf0
            - 2 * self._equation.Ut_start_x(self._tau * (k + 1)) * self._h
        ) / denominator_l

        self._answer[k + 1][-1] = (
            4 * self._answer[k + 1][-2] * self._equation.alfn
            - self._answer[k + 1][-3] * self._equation.alfn
            + 2 * self._equation.Ut_finish_x(self._tau * (k + 1)) * self._h
        ) / denominator_r

    def _two_points_second_order(self, k):
        denominator_l = (
            2 * self._equation.alf0 * self._tau**2 * self._equation.a
            - self._equation.alf0 * self._tau**2 * self._equation.c * self._h**2
            + self._equation.alf0 * self._tau * self._equation.e * self._h**2
            + self._equation.alf0 * self._h**2
            - 2 * self._equation.bet0 * self._tau**2 * self._equation.a * self._h
            + self._equation.bet0 * self._tau**2 * self._equation.b * self._h**2
        )
        denominator_r = (
            2 * self._equation.alfn * self._tau**2 * self._equation.a
            - self._equation.alfn * self._tau**2 * self._equation.c * self._h**2
            + self._equation.alfn * self._tau * self._equation.e * self._h**2
            + self._equation.alfn * self._h**2
            + 2 * self._equation.betn * self._tau**2 * self._equation.a * self._h
            + self._equation.betn * self._tau**2 * self._equation.b * self._h**2
        )

        self._answer[k + 1][0] = (
            (
                self._equation.alf0 * self._tau * self._equation.e * self._h**2
                + 2 * self._equation.alf0 * self._h**2
            )
            * self._answer[k][0]
            - self._equation.alf0 * self._h**2 * self._answer[k - 1][0]
            + 2
            * self._equation.alf0
            * self._tau**2
            * self._equation.a
            * self._answer[k + 1][1]
            + self._equation.alf0
            * self._tau**2
            * self._h**2
            * self.equation.fi(self._x[0], self._tau * (k + 1))
            + (
                self._tau**2 * self._equation.b * self._h**2
                - 2 * self._tau**2 * self._h * self._equation.a
            )
            * self._equation.Ut_start_x(self._tau * (k + 1))
        ) / denominator_l

        self._answer[k + 1][-1] = (
            (
                self._equation.alfn * self._tau * self._equation.e * self._h**2
                + 2 * self._equation.alfn * self._h**2
            )
            * self._answer[k][-1]
            - self._equation.alfn * self._h**2 * self._answer[k - 1][-1]
            + 2
            * self._equation.alfn
            * self._tau**2
            * self._equation.a
            * self._answer[k + 1][-2]
            + self._equation.alfn
            * self._tau**2
            * self._h**2
            * self.equation.fi(self._x[-1], self._tau * (k + 1))
            + (
                self._tau**2 * self._equation.b * self._h**2
                + 2 * self._tau**2 * self._h * self._equation.a
            )
            * self._equation.Ut_finish_x(self._tau * (k + 1))
        ) / denominator_r

    def _calc_answer(self):
        # print(self._equation.a*self._tau**2 / self._h**2)
        # print("calc")
        self._answer = np.zeros((self._resolution_t, self._resolution_x))
        self._answer[0, :] = [self._equation.Ux_start_t(elem) for elem in self._x]
        if self._k1_method == 1:
            self._k1_method_two_point()
        elif self._k1_method == 2:
            self._k1_method_taylor_series()
        for k in range(1, self._resolution_t - 1):
            self._calc_midle_points(k)

            if self._lr_method == 1:
                self._two_points_first_order(k)
            elif self._lr_method == 2:
                self._three_points_second_order(k)
            elif self._lr_method == 3:
                self._two_points_second_order(k)
