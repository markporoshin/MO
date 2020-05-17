import numpy as np
from enum import Enum


def f(x):
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3]


def df(x):
    return np.array([2 * x[0], 2 * x[1], 2 * x[2], 2 * x[3]])


C = np.array([[0.1, 0, 0, 0], [0, 0.1, 0, 0]])
d = np.array([0, 0])

F = np.array([[0.01, 0.01, 0.01, 0.01], [-0.01, -0.01, -0.01, -0.01]])
g = np.array([1, 1])


class State(Enum):
    PROCESS = 1
    ANSWER = 2


def find_step(f, a0, l, F2, g2, sk, xk):
    a_more = a0
    a_less = a0
    while True:
        xk_l = xk + a_less * sk
        xk_m = xk + a_more * sk
        if f(xk_l) < f(xk) and (F2.dot(xk_l) <= g2).all():
            return a_less
        if f(xk_m) < f(xk) and (F2.dot(xk_m) <= g2).all():
            return a_more
        a_less = a_less / l
        a_more = a_more * l


def find_dir(df, C, d, F1, f1_indexes, g, xk):
    if F1.shape[0] != 0 and F1.shape[1] != 0:
        a_mtr = np.concatenate(C, F1)
    else:
        a_mtr = C.copy()
    E = np.identity(F1.shape[1])
    if a_mtr.shape == (0, 0):
        pk_mtr = E
    else:
        pk_mtr = E - a_mtr.transpose().dot(np.linalg.inv(a_mtr.dot(a_mtr.transpose()))).dot(a_mtr)
    sk = -pk_mtr.dot(df(xk))
    if (sk != 0).any():
        return sk, State.PROCESS
    elif a_mtr.shape == (0, 0) and (sk == 0).all():
        return xk, State.ANSWER
    else:
        w = - np.linalg.inv(a_mtr * a_mtr.transpose()).dot(a_mtr).dot(df(xk))
        u = w[f1_indexes:]
        if (u > 0).all():
            return xk, State.ANSWER
        else:
            u_j = np.where(u == np.amin(u))
            F1 = np.delete(F1, u_j, axis=0)
            f1_indexes.remove(f1_indexes.index(u_j))
            return find_dir(df, C, d, F1, f1_indexes, g, xk)


def method(f, df, C, d, F, g, xk, a0, l):
    state = State.PROCESS
    while state != State.ANSWER:
        # prepare F1
        f1_indexes = []
        for ind in range(C.shape[0]):
            if F[ind].dot(xk) == d[ind]:
                f1_indexes.append(ind)
        F1 = F[f1_indexes]

        # find direction
        sk, state = find_dir(df, C, d, F1, f1_indexes, g, xk)
        if state == State.ANSWER:
            return xk

        # find step
        f2_indexes = [_ for _ in range(C.shape[0]) if _ not in f1_indexes]
        F2 = F[f2_indexes]
        g2 = g[f2_indexes]
        ak = find_step(f, a0, l, F2, g2, sk, xk)

        # go
        xk = xk + ak * sk


xk = method(f, df, C, d, F, g, np.array([3, 3, 3, 3]), 1, 0.7)
print(xk)