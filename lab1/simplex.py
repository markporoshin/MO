import numpy as np
from sv import SV
import random

def list_diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif


def count_tetta(x, u, N):
    ik = 0
    for i in N:
        if u[i] > 0:
            ik = i
            break
    tetta = x[ik] / u[ik]
    for i in N:
        if u[i] > 0 and tetta > x[i] / u[i]:
            tetta = x[i] / u[i]
            ik = i
    return tetta, ik


# we don't need j, it in u
def next_b(u, B, i, N):
    F = np.identity(len(N)) # mb ERROR
    for row in range(0, len(N)):
        F[row][N.index(i)] = -u[N[row]] / u[i]
    F[N.index(i)][N.index(i)] = 1 / u[i]
    return F.dot(B)


def change_basis(sv, A):
    # Nk = sv.N.copy()
    i = random.choice(list_diff(sv.N, sv.N_plus))
    j = random.choice(sv.L)

    # sv.L.pop(sv.L.index(j))
    # sv.L.append(i)
    # sv.N.pop(sv.N.index(i))
    # sv.N.append(j)

    sv.L[sv.L.index(j)] = i
    sv.N[sv.N.index(i)] = j
    sv.A[:, sv.N.index(j)] = A[:, j]

    # sv.B = next_b(A[:, j], sv.B, sv.N.index(j), Nk)
    # sv.B = np.linalg.inv(sv.A)
    while np.linalg.det(sv.A) == 0:
        i = random.choice(list_diff(sv.N, sv.N_plus))
        j = random.choice(sv.L)
        # sv.L.pop(sv.L.index(j))
        # sv.L.append(i)
        # sv.N.pop(sv.N.index(i))
        # sv.N.append(j)

        sv.L[sv.L.index(j)] = i
        sv.N[sv.N.index(i)] = j
        sv.A[:, sv.N.index(j)] = A[:, j]
        # sv.B = next_b(A[:, j], sv.B, i, Nk)
    sv.B = np.linalg.inv(sv.A)

def simplex(A, b, c, sv, M, N):
    i = 1
    # print("------step 1------")
    state = 'run'
    # print(sv.x)
    # print("target function:" + str(c.dot(sv.x)))
    # print("-------------------------")
    while state == 'run' or state == 'next_sv':
        # print("------step " + str(i)+"-------")
        state = step(A, b, c, sv, M, N)
        print("next sv" + str(sv.x))
        print("target function:" + str(c.dot(sv.x)))
        # print("-------------------------")
        i += 1
    return state


def build_next_x(sv, jk, u, A, Nk, N):
    tetta, ik = count_tetta(sv.x, u, Nk)
    sv.set_x(sv.x - tetta * u, N)
    sv.A[:, Nk.index(ik)] = A[:, jk]
    sv.N[sv.N.index(ik)] = jk
    sv.L[sv.L.index(jk)] = ik
    sv.B = next_b(u, sv.B, ik, Nk)
    # print("Bk+1:")
    # print(sv.B)
    #
    #
    # print("ik: "+ str(ik) + " jk " + str(jk))
    # print(sv.N)
    # print(sv.A)
    if not (sv.B == np.linalg.inv(sv.A)).all():
        # print("inverse matrix problems")
        # print(sv.A * sv.B)
        sv.B = np.linalg.inv(sv.A)


def step(A, b, c, sv, M, N):
    # check section
    if (A[:, sv.N] != sv.A).any():
        print("ERROR: ...")
    # if (A.dot(sv.x) != b).any():
    #     print("ERROR: incorrect support vector")
    #     print(A.dot(sv.x))
    # if (sv.A.dot(sv.B) != np.identity(len(M))).any():
    #     print("ERROR: incorrect inverse matrix")
    #     sv.B = np.linalg.inv(sv.A)

    Nk = sv.N.copy()
    Lk = sv.L.copy()
    Bk = sv.B.copy()
    y = sv.B.transpose().dot(c[Nk])
    d_L = c[Lk] - A[:, Lk].transpose().dot(y.transpose())
    d_N = c[Nk] - A[:, Nk].transpose().dot(y.transpose())
    if (d_N != 0).all():
        print("ERROR in d[Nk]")
    d = np.zeros(len(N))
    d[Lk] = d_L
    d[Nk] = d_N
    # print("d:")
    # print(d)
    if (d_L >= 0).all():
        return "optimal"

    jk = 0
    for e in N:
        if d_L[e] < 0:
            jk = Lk[e]
            break

    u = np.zeros(len(N))
    u[Nk] = Bk.dot(A[M, jk])
    u[Lk] = 0
    u[jk] = -1
    # print("u:")
    # print(u)
    if (u[Nk] <= 0).all():
        x = sv.x - 1000 * u
        if (A.dot(x) != b).any():
            print("ERROR: incorrect support vector")
        return "unbounded"

    if sorted(sv.N) == sorted(sv.N_plus):
        build_next_x(sv, jk, u, A, Nk, N)
        return "next_sv"
    else:
        if (u[list_diff(Nk, sv.N_plus)] <= 0).all():
            # print("next sv")
            # print(sv.x)
            build_next_x(sv, jk, u, A, Nk, N)
            return "next_sv"
        else:
            change_basis(sv, A)
            # print(sv.N)
            # print(sv.L)
            # print(sv.A)
            return step(A, b, c, sv, M, N)




