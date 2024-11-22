import numpy as np
import solver


def trid_alg(matrix, vec):
    s = len(vec)
    ans = np.zeros(s)
    p = np.zeros(s)
    q = np.zeros(s)

    p[0] = -matrix[0][1] / matrix[0][0]
    q[0] = vec[0] / matrix[0][0]

    for i in range(1, s - 1):
        p[i] = -matrix[i][i + 1] / (matrix[i][i] + matrix[i][i - 1] * p[i - 1])
        q[i] = (vec[i] - matrix[i][i - 1] * q[i - 1]) / (
            matrix[i][i] + matrix[i][i - 1] * p[i - 1]
        )

    p[s - 1] = 0
    q[s - 1] = (vec[s - 1] - matrix[s - 1][s - 2] * q[s - 2]) / (
        matrix[s - 1][s - 1] + matrix[s - 1][s - 2] * p[s - 2]
    )

    ans[s - 1] = q[s - 1]

    for i in range(s - 2, -1, -1):
        ans[i] = p[i] * ans[i + 1] + q[i]

    return ans


class Hyperbolic_implicit_solver(solver.Solver):
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
        fraction_a = self._tau**2 * self._equation.a / self._h**2
        fraction_b = self._tau**2 * self._equation.b / (4 * self._h)
        fraction_c = self._tau**2 * self._equation.c / 2
        fraction_e = self._tau**2 * self._equation.e / 2

        for i in range(1, self._resolution_x - 1):
            self._answer[1][i] = (
                (fraction_c - fraction_a + 1) * self._answer[0][i]
                + (self._tau - fraction_e) * self._equation.Ux_start_t_dt(self._x[i])
                + (fraction_a / 2 + fraction_b) * self._answer[0][i + 1]
                + (fraction_a / 2 - fraction_b) * self._answer[0][i - 1]
                + self.equation.fi(self._x[i], 0) * self._tau**2 / 2
            )

    def _set_mean_equations_system(self, k):
        denominator = (
            2 * self._tau**2 * self._equation.a
            - self._tau**2 * self._equation.c * self._h**2
            + self._tau * self._equation.e * self._h**2
            + self._h**2
        )
        fraction_a = self._tau**2 * self._equation.a
        fraction_b = self._tau**2 * self._equation.b * self._h
        fraction_e = self._tau * self._equation.e * self._h**2
        for i in range(1, self._resolution_x - 1):
            self._new_line_mat[i][i] = denominator
            self._new_line_mat[i][i + 1] = -(fraction_a + fraction_b / 2)
            self._new_line_mat[i][i - 1] = -(fraction_a - fraction_b / 2)
            self._new_line_vec[i] = (
                (fraction_e + 2 * self._h**2) * self._answer[k][i]
                - self._h**2 * self._answer[k - 1][i]
                + self._h**2
                * self._tau**2
                * self.equation.fi(self._x[i], self._tau * (k + 1))
            )

    def _two_points_first_order(self, k):
        denominator_l = self._equation.alf0 - self._equation.bet0 * self._h
        denominator_r = self._equation.betn * self._h + self._equation.alfn

        self._new_line_mat[0][0] = denominator_l
        self._new_line_mat[0][1] = -self._equation.alf0
        self._new_line_vec[0] = (
            -self._equation.Ut_start_x(self._tau * (k + 1)) * self._h
        )

        self._new_line_mat[-1][-1] = denominator_r
        self._new_line_mat[-1][-2] = -self._equation.alfn
        self._new_line_vec[-1] = (
            self._equation.Ut_finish_x(self._tau * (k + 1)) * self._h
        )

    def _three_points_second_order(self, k):
        denominator_l = 3 * self._equation.alf0 - 2 * self._equation.bet0 * self._h
        denominator_r = 3 * self._equation.alfn + 2 * self._equation.betn * self._h

        self._new_line_mat[0][0] = denominator_l
        self._new_line_mat[0][1] = -4 * self._equation.alf0
        self._new_line_mat[0][2] = self._equation.alf0
        self._new_line_vec[0] = (
            -2 * self._equation.Ut_start_x(self._tau * (k + 1)) * self._h
        )

        self._new_line_mat[-1][-1] = denominator_r
        self._new_line_mat[-1][-2] = -4 * self._equation.alfn
        self._new_line_mat[-1][-3] = self._equation.alfn
        self._new_line_vec[-1] = (
            2 * self._equation.Ut_finish_x(self._tau * (k + 1)) * self._h
        )

        coef_0 = self._new_line_mat[0][2] / self._new_line_mat[1][2]
        coef_n = self._new_line_mat[-1][-3] / self._new_line_mat[-2][-3]

        self._new_line_mat[0] -= self._new_line_mat[1] * coef_0
        self._new_line_vec[0] -= self._new_line_vec[1] * coef_0
        self._new_line_mat[-1] -= self._new_line_mat[-2] * coef_n
        self._new_line_vec[-1] -= self._new_line_vec[-2] * coef_n

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

        self._new_line_mat[0][0] = denominator_l
        self._new_line_mat[0][1] = (
            -2 * self._equation.alf0 * self._tau**2 * self._equation.a
        )
        self._new_line_vec[0] = (
            (
                self._equation.alf0 * self._tau * self._equation.e * self._h**2
                + 2 * self._equation.alf0 * self._h**2
            )
            * self._answer[k][0]
            - self._equation.alf0 * self._h**2 * self._answer[k - 1][0]
            + self._equation.alf0
            * self._tau**2
            * self._h**2
            * self.equation.fi(self._x[0], self._tau * (k + 1))
            + (
                self._tau**2 * self._equation.b * self._h**2
                - 2 * self._tau**2 * self._h * self._equation.a
            )
            * self._equation.Ut_start_x(self._tau * (k + 1))
        )

        self._new_line_mat[-1][-1] = denominator_r
        self._new_line_mat[-1][-2] = (
            -2 * self._equation.alfn * self._tau**2 * self._equation.a
        )
        self._new_line_vec[-1] = (
            (
                self._equation.alfn * self._tau * self._equation.e * self._h**2
                + 2 * self._equation.alfn * self._h**2
            )
            * self._answer[k][-1]
            - self._equation.alfn * self._h**2 * self._answer[k - 1][-1]
            + self._equation.alfn
            * self._tau**2
            * self._h**2
            * self.equation.fi(self._x[-1], self._tau * (k + 1))
            + (
                self._tau**2 * self._equation.b * self._h**2
                + 2 * self._tau**2 * self._h * self._equation.a
            )
            * self._equation.Ut_finish_x(self._tau * (k + 1))
        )

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
            self._new_line_mat = np.zeros((self._resolution_x, self._resolution_x))
            self._new_line_vec = np.zeros(self._resolution_x)
            self._set_mean_equations_system(k)
            if self._lr_method == 1:
                self._two_points_first_order(k)
            elif self._lr_method == 2:
                self._three_points_second_order(k)
            elif self._lr_method == 3:
                self._two_points_second_order(k)
            new_line = trid_alg(self._new_line_mat, self._new_line_vec)
            self._answer[k + 1] = new_line
