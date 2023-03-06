import numpy as np
import matplotlib.pyplot as plt

class Laguerre:

    def __init__(self, t, n, beta, sigma, fun, a, b, number_of_points, eps, T):
        self.t = t
        self.n = n
        self.beta = beta
        self.sigma = sigma
        self.fun = fun
        self.a = a
        self.b = b
        self.number_of_points = number_of_points
        self.eps = eps
        self.T = T

    def laguerre_func(self, t, n, beta=2, sigma=4):

        laguerre_result_0 = np.sqrt(sigma) * (np.exp(-(beta * t) / 2))
        laguerre_result_1 = laguerre_result_0 * (1 - (sigma * t))

        if n == 0:
            return laguerre_result_0
        elif n == 1:
            return laguerre_result_1
        else:
            for i in range(2, n + 1):
                laguerre_result_n = ((2 * i - 1 - t * sigma) / i) * laguerre_result_1 - (
                            (i - 1) / i) * laguerre_result_0
                laguerre_result_0 = laguerre_result_1
                laguerre_result_1 = laguerre_result_n
            return laguerre_result_n

    def laguerre_tabulation(self, t, n, beta, sigma):
        laguerre_result_0 = self.laguerre_func(self.t, 0, self.beta, self.sigma)
        laguerre_result_1 = self.laguerre_func(self.t, 1, self.beta, self.sigma)
        list_of_lag = [laguerre_result_0, laguerre_result_1]
        list_of_t_values = [i for i in range(0, self.t + 1)]
        for i in range(2, t + 1):
            laguerre_result_n = self.laguerre_func(self.t, i, self.beta, self.sigma)
            laguerre_result_0 = laguerre_result_1
            laguerre_result_1 = laguerre_result_n
            list_of_lag.append(laguerre_result_n)
        return list_of_lag, list_of_t_values

    def function_graph(self, fun, a, b, number_of_points, n, sigma, beta):
        fig = plt.figure(figsize=(7, 5))
        ax = fig.gca()

        x = np.linspace(self.a, self.b, self.number_of_points)

        for self.n in range(21):
            y = np.array([func_val for func_val in fun(self.x, self.n, self.beta, self.sigma)])
            ax.plot(x, y)

        ax.axhline(color='grey')
        ax.axvline(color='grey')
        ax.set_xlim(self.a, self.b)
        ax.set_xlabel('x')
        ax.set_ylabel('lag_func(x)')
        ax.grid()
        plt.show()

    def rectangle_integral(self, fun, a, T, number_of_points, eps, N, beta, sigma):
        alpha = self.sigma - self.beta
        delta = abs((self.T - self.a) / (self.number_of_points - 1))
        half_delta = delta / 2
        points = np.linspace(half_delta, self.T - half_delta, self.number_of_points - 1)
        arr1 = np.array([np.sum(
            [fun(i) * self.laguerre_func(i, self.n, self.beta, self.sigma) * (np.e ** ((-1) * alpha * i)) * delta for i
             in points]) for self.n in range(0, self.N + 1)])
        self.number_of_points *= 2
        delta = abs((self.T - self.a) / (self.number_of_points - 1))
        half_delta = delta / 2
        points = np.linspace(half_delta, self.T - half_delta, self.number_of_points - 1)
        arr2 = np.array([np.sum(
            [fun(i) * self.laguerre_func(i, self.n, self.beta, self.sigma) * (np.e ** ((-1) * alpha * i)) * delta for i
             in points]) for self.n in range(0, self.N + 1)])
        while (abs(arr2 - arr1) > eps).any():
            arr1 = arr2
            self.number_of_points *= 2
            delta = abs((self.T - self.a) / (self.number_of_points - 1))
            half_delta = delta / 2
            points = np.linspace(half_delta, self.T - half_delta, self.number_of_points - 1)
            arr2 = np.array([np.sum(
                [fun(i) * self.laguerre_func(i, self.n, self.beta, self.sigma) * (np.e * ((-1) * alpha * i)) * delta for
                 i in points]) for self.n in range(0, self.N + 1)])
        return arr2

    def seq_sum(self, N, beta, sigma, t):
        sm = 0
        seq = np.array(self.N)
        for n in range(len(self.N)):
            if seq[n] != 0:
                sm += seq[n] * self.laguerre_func(self.t, self.n, self.beta, self.sigma)
            else:
                break
        return sm


def f(x):
    if (x < 2 * np.pi) and (x >= 0):
        return np.sin(x - np.pi / 2) + 1
    elif (x >= 2 * np.pi):
        return 0

t, n, beta, sigma, fun, a, b, number_of_points, eps, T = 5, 5, 2, 4, f, 5, 29, 100, 0.01, 15
lag = Laguerre(t, n, beta, sigma, fun, a, b, number_of_points, eps, T)
lag.laguerre_tabulation(t, n, beta, sigma)