import numpy as np
from math import sin, cos


def f(x):
    a = 5
    b = 4
    c = 5
    return x[0] * x[0] + a * x[1] * x[1] + sin(b * x[0] + c * x[1]) + 3 * x[0] + 2 * x[1]


def df(x):
    a = 5
    b = 4
    c = 5
    return np.array([2 * x[0] + b * cos(b * x[0] + c * x[1]) + 3,
                     2 * a * x[1] + c * cos(b * x[0] + c * x[1]) + 2])


def grad(df, xk, l, e):
    xk_next = xk - l * df(xk)
    while np.linalg.norm(xk_next - xk) > e:
        xk = xk_next
        xk_next = xk - l * df(xk)
    return xk_next

print(grad(df, np.array([0, 0]), 1e-2, 1e-2))
print(f(grad(df, np.array([0, 0]), 1e-2, 1e-2)))