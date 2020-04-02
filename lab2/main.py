import numpy as np

f = open('size.txt', 'r')
line1 = f.read()
n = int(line1.split(' ')[0])
m = int(line1.split(' ')[1])

C = np.loadtxt('C.txt')
a = np.loadtxt('a.txt')
b = np.loadtxt('b.txt')

# checked
def northwest_angle(a, b, n, m):
    a_copy = a.copy()
    b_copy = b.copy()

    # current[0] - row
    # current[1] - col
    # current[2] - value
    current = [0, 0, 0]
    x = np.zeros((n, m))
    x[:, :] = -1
    while current[0] < n and current[1] < m:
        v = min(a_copy[current[0]], b_copy[current[1]])
        a_copy[current[0]] -= v
        b_copy[current[1]] -= v
        current[2] = v
        x[current[0]][current[1]] = v
        if a_copy[current[0]] > 0:
            current[1] += 1
        else:
            current[0] += 1
    return x


def count_potential(C, x, n, m):
    A = np.zeros((n + m, n + m))
    b = np.zeros(n + m)
    row = 1
    A[0][0] = 1
    b[0] = 0

    for r in range(n):
        for c in range(m):
            if x[r][c] != -1:
                A[row][r] = 1
                A[row][n+c] = 1
                b[row] = C[r][c]
                row += 1
    u_v = np.linalg.solve(A, b)
    u = u_v[0: n]
    v = u_v[n: n + m]
    return u, v


def count_deltas(u, v, x, C, n, m):
    deltas = np.zeros(n, m)
    for r in range(n):
        for c in range(m):
            if x[r][c] == 0:
                deltas[r][c] = C[r][c] - (u[r] + v[c])
    return deltas


def find_row(x, point, stack):
    pass


def find_col(x, point, stack):
    pass


def find_cycle(x, point):
    stack = []

    find_col(x, point, stack)

    return stack


def improve_plan(x, u, v, a, b, C, deltas):
    point = deltas.argmin()
    cycle = find_cycle(x, point)

    pass


def method(a, b, C, n, m):
    x = northwest_angle(a, b, n, m)
    u, v = count_potential(C, x, n, m)
    deltas = count_deltas(u, v, x, C, n, m)
    while (deltas < 0).any():
        improve_plan(x, u, v, a, b, C, deltas)
        u, v = count_potential(C, x, n, m)
        deltas = count_deltas(u, v, x, C, n, m)


method(a, b, C, n, m)
