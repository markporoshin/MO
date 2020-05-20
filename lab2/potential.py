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
    deltas = np.zeros((n, m))
    for r in range(n):
        for c in range(m):
            if x[r][c] == -1:
                deltas[r][c] = C[r][c] - (u[r] + v[c])
    return deltas


def find_row(x, point, stack, n, m, src):
    if src == point and len(stack) > 0:
        return True
    row = point[0]
    for c in range(m):
        if (x[row][c] != -1 or (row, c) == src) and (row, c) != point:
            if find_col(x, (row, c), stack, n, m, src):
                stack.append(point)
                return True
    return False


def find_col(x, point, stack, n, m, src, flag=True):
    if src == point and flag:
        return True
    col = point[1]
    for r in range(n):
        if (x[r][col] != -1 or (r, col) == src) and (r, col) != point:
            if find_row(x, (r, col), stack, n, m, src):
                stack.append(point)
                return True
    return False

def find_cycle(x, point, n, m):
    stack = []
    find_col(x, point, stack, n, m, point, False)
    return stack


def recount_x(x, cycle):
    i = 0
    minus = []
    plus = []
    for el in cycle:
        if i % 2 == 0:
            minus.append(el)
        else:
            plus.append(el)
        i += 1
    minimum = min([x[el[0]][el[1]] for el in minus])
    flag = True
    for el in minus:
        x[el[0]][el[1]] -= minimum
        if x[el[0]][el[1]] == 0 and flag:
            x[el[0]][el[1]] = -1
            flag = False
    for el in plus:
        if x[el[0]][el[1]] == -1:
            x[el[0]][el[1]] += 1 + minimum
        else:
            x[el[0]][el[1]] += minimum


def improve_plan(x, a, b, C, deltas, n, m):
    point = np.unravel_index(np.argmin(deltas, axis=None), deltas.shape)
    cycle = find_cycle(x, point, n, m)
    recount_x(x, cycle)
    pass


def method(a, b, C, n, m):
    x = northwest_angle(a, b, n, m)
    print(x)
    u, v = count_potential(C, x, n, m)
    deltas = count_deltas(u, v, x, C, n, m)
    while (deltas < 0).any():
        print("step")
        improve_plan(x, a, b, C, deltas, n, m)
        u, v = count_potential(C, x, n, m)
        deltas = count_deltas(u, v, x, C, n, m)
    return x


x = method(a, b, C, n, m)
print()
print(x)
