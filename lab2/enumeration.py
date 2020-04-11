import numpy as np

def next_set(set, n, m):
    k = m
    for i in reversed(range(k)):
        if set[i] < n-k+i+1:
            set[i] += 1
            for j in range(i+1, k):
                set[j] = set[j-1] + 1
            return True
    return False


def count_sv(inds, A, b, m):
    if np.linalg.det(A[:, inds]) == 0:
        return []

    sv = np.zeros(m)
    sv[inds] = np.linalg.solve(A[:,inds], b)
    if (sv < 0).any():
        return []
    return sv


f = open('../lab1/input.txt', 'r')
line1 = f.read()
n = int(line1.split(' ')[0])
m = int(line1.split(' ')[1])

A = np.loadtxt('../lab1/A.txt')
c = np.loadtxt('../lab1/c.txt')
b = np.loadtxt('../lab1/b.txt')


set = [_ + 1 for _ in range(m)]
min = set[0:n]
cur_sv = count_sv([i-1 for i in min], A, b, m)
target_min = float('inf')
if cur_sv != []:
    target_min = c.dot(cur_sv)


while next_set(set, m, n):
    indexes = set[0:n]
    if indexes == [1, 3, 4, 6, 12, 16, 17, 20]:
        print("equal")

    cur_sv = count_sv([i-1 for i in indexes], A, b, m)
    if cur_sv == []:
        continue
    target = c.dot(cur_sv)
    if target_min > target:
        min = indexes
        target_min = target
        print(cur_sv)
        print(target)
print("min value %s" % target_min)
print(count_sv([i-1 for i in min], A, b, m))
