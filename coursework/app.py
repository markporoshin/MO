import numpy as np
from math import sqrt, fabs


def golden_ratio(fun, a, b, E):
    print("golden_section_method")
    K = (1 + sqrt(5)) / 2
    lam = b - (b - a) / K
    mu = a + (b - a) / K
    y1 = fun(lam)
    y2 = fun(mu)
    counter = 2
    n = 0
    while fabs(b - a) > E:
        print("step %s: lam-%s mu-%s" % (n, lam, mu))
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
            lam = lam = b - (b - a) / K
            y1 = fun(lam)
            counter += 1
        n += 1
    return (b + a) / 2, counter


def uniform_search(fun, a, b, E, n=5, counter=0):
    if b - a < E:
        return (b - a) / 2, counter
    step = (b - a) / n
    xn = a + step
    min_x = a
    min_f = fun(a)
    while xn < b:
        counter += 1
        f = fun(xn)
        if f < min_f:
            min_x = xn
            min_f = f
        xn += step
    counter += 2
    if fun(min_x - step) < fun(min_x + step):
        return uniform_search(fun, min_x - step, min_x, E, counter)
    return uniform_search(fun, min_x, min_x + step, E, counter)


def steepest_descent_method(fun, dfun, xk, E, minimize):
    s0 = dfun(xk)
    phi = lambda a: fun(xk + s0 * a)
    xk_next = xk - minimize(phi, 0, 1, E) * s0

    while np.linalg.norm(xk_next - xk) > E:
        sk = dfun(xk)
        phi = lambda a: fun(xk + sk * a)
        xk_next = xk - minimize(phi, 0, 1, E) * s0

    return xk_next, fun(xk_next)










