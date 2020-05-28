import numpy as np
from math import sqrt, fabs, sin, cos
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches


def golden_ratio(fun, a, b, E):
    """
    find a minimum of one-dimensional and unimodal function using golden ratio method
    :param fun: the function to be minimized
    :param a: left border
    :param b: right border
    :param E: accuracy
    :return: x - point of minimum, f_x - value of function in x, counter - number of function call
    """
    K = (1 + sqrt(5)) / 2
    lam = b - (b - a) / K
    mu = a + (b - a) / K
    y1 = fun(lam)
    y2 = fun(mu)
    counter = 2
    n = 0
    while fabs(b - a) > E:
        if y1 > y2:
            a = lam
            lam = mu
            y1 = y2
            mu = a + (b - a) / K
            y2 = fun(mu)
            counter += 1
        else:
            b = mu
            mu = lam
            y2 = y1
            lam = b - (b - a) / K
            y1 = fun(lam)
            counter += 1
        n += 1
    return (b + a) / 2, counter


def uniform_search(fun, a, b, E, n=3, counter=0):
    """
    find a minimum of one-dimensional and unimodal function using uniform search method
    :param fun: the function to be minimized
    :param a: left border
    :param b: right border
    :param E: accuracy
    :return: x - point of minimum, f_x - value of function in x, counter - number of function call
    """
    if b - a < E:
        return (b + a) / 2, counter
    step = (b - a) / n
    xn = a + step
    min_x = a
    min_f = fun(a)
    while xn <= b:
        counter += 1
        f = fun(xn)
        if f < min_f:
            min_x = xn
            min_f = f
        xn += step
    counter += 2
    if fun(min_x - step) < fun(min_x + step):
        return uniform_search(fun, min_x - step, min_x, E, n, counter)
    return uniform_search(fun, min_x, min_x + step, E, n, counter)


def steepest_descent_method(fun, dfun, xk, E, minimize):
    """
        find a minimum function using uniform search steepest descent method
        :param fun: the function to be minimized
        :type fun: function(np.array)
        :param dfun: the differential of function to be minimized
        :type dfun: function(np.array)
        :param xk: start point
        :type xk: np.array
        :param minimize:
        :type minimize: function(function(np.array), double, double, double)
        :param E: accuracy
        :return: x - point of minimum, f_x - value of function in x, counter - number of function call
        """
    counter = 0
    s0 = -dfun(xk)
    p = lambda a: fun(xk + s0 * a)
    a0, f_appeals = minimize(p, 0, 1, E)
    counter += f_appeals
    xk_next = xk + a0 * s0

    while np.linalg.norm(xk_next - xk) > E:
        xk = xk_next
        sk = -dfun(xk)
        p = lambda a: fun(xk + sk * a)
        ak, f_appeals = minimize(p, 0, 1, E)
        counter += f_appeals
        xk_next = xk + ak * sk

    return xk_next, fun(xk_next), counter


def f1(x):
    return x[0] * x[0] + x[1] * x[1] + sin(x[0] * x[1])


def df1(x):
    return np.array([2 * x[0] + x[1] * cos(x[0] * x[1]), 2 * x[1] + x[0] * cos(x[0] * x[1])])


def f2(x):
    a = 5
    b = 4
    c = 5
    return x[0] * x[0] + a * x[1] * x[1] + sin(b * x[0] + c * x[1]) + 3 * x[0] + 2 * x[1]


def df2(x):
    a = 5
    b = 4
    c = 5
    return np.array([2 * x[0] + b * cos(b * x[0] + c * x[1]) + 3,
                     2 * a * x[1] + c * cos(b * x[0] + c * x[1]) + 2])


def f3(x):
    return 3 * x[0] * x[0] + 7 * x[1] * x[1] + x[2] * x[2] + (x[3] - 9) * (x[3] - 9) + x[0] * x[1] * sin(x[0] + x[1])


def df3(x):
    return np.array([2 * x[0] + x[0] * x[1] * cos(x[0] + x[1]) + x[1] * sin(x[0] + x[1]),
                     2 * x[1] + x[0] * x[1] * cos(x[0] + x[1]) + x[0] * sin(x[0] + x[1]),
                     2 * x[2],
                     2 * (x[3] - 9)])


x03 = np.array([6, 1, -1, 3])

us3 = lambda fun, a, b, E: uniform_search(fun, a, b, E, n=3)
us2 = lambda fun, a, b, E: uniform_search(fun, a, b, E, n=2)
us10 = lambda fun, a, b, E: uniform_search(fun, a, b, E, n=10)

E = 1
for i in range(8):
    _, _, numOfCallsGR = steepest_descent_method(f3, df3, x03, E, golden_ratio)
    _, _, numOfCallsUS3 = steepest_descent_method(f3, df3, x03, E, us3)
    _, _, numOfCallsUS2 = steepest_descent_method(f3, df3, x03, E, us2)
    _, _, numOfCallsUS10 = steepest_descent_method(f3, df3, x03, E, us10)
    plot.plot(i, numOfCallsGR, 'r.')
    plot.plot(i, numOfCallsUS3, 'g.')
    plot.plot(i, numOfCallsUS2, 'b.')
    plot.plot(i, numOfCallsUS10, 'y.')
    E = E / 10
plot.xlabel('i')
plot.ylabel('calls')
legends = []
legends.append(mpatches.Patch(color='r', label='золотое сечение'))
legends.append(mpatches.Patch(color='g', label='равномерный поиск n=3'))
legends.append(mpatches.Patch(color='b', label='равномерный поиск n=2'))
legends.append(mpatches.Patch(color='y', label='равномерный поиск n=10'))
plot.legend(handles=legends)
plot.show()
# print(f"{x}, {f_x}, {numOfCalls}")
