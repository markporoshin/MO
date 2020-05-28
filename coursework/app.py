import numpy as np
from math import sqrt, fabs, sin, cos


def golden_ratio(fun, a, b, E):
    # print("golden_section_method")
    K = (1 + sqrt(5)) / 2
    lam = b - (b - a) / K
    mu = a + (b - a) / K
    y1 = fun(lam)
    y2 = fun(mu)
    # print(f"[{lam}, {mu}]")
    counter = 2
    n = 0
    while fabs(b - a) > E:
        # print("step %s: lam-%s mu-%s" % (n, lam, mu))
        if y1 > y2:
            a = lam
            lam = mu
            y1 = y2
            mu = a + (b - a) / K
            y2 = fun(mu)
            counter += 1
            # print(f"[{lam}, {mu}]")
        else:
            b = mu
            mu = lam
            y2 = y1
            lam = b - (b - a) / K
            y1 = fun(lam)
            counter += 1
            # print(f"[{lam}, {mu}]")
        n += 1
    return (b + a) / 2, counter


def uniform_search(fun, a, b, E, n=2, counter=0):
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


def phi(fun, x, a, s):
    return fun(x + s * a)


def steepest_descent_method(fun, dfun, xk, E, minimize):
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


def f(x):
    return x[0] * x[0] + x[1] * x[1] + sin(x[0] * x[1])


def df(x):
    return np.array([2 * x[0] + x[1] * cos(x[0] * x[1]), 2 * x[1] + x[0] * cos(x[0] * x[1])])


x, f_x, numOfCalls = steepest_descent_method(f, df, np.array([2, 2]), 1e-2, golden_ratio)
print(f"{x}, {f_x}, {numOfCalls}")









