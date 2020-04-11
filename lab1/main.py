import numpy as np
from lab1.sv import SV
from lab1.simplex import simplex


def list_diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif


f = open('input.txt', 'r')
line1 = f.read()
ah = int(line1.split(' ')[0])
aw = int(line1.split(' ')[1])
A = np.loadtxt('A.txt')
b = np.loadtxt('b.txt')
c = np.loadtxt('c.txt')
M = [i for i in range(ah)]
N = [i for i in range(aw)]

s_A = np.append(A, np.identity(ah), axis=1)
s_x = np.append(np.zeros(shape=(aw, 1)), b)
s_c = np.append(np.zeros(shape=(aw, 1)), np.ones(shape=(ah, 1)))
s_ah = ah
s_aw = aw + ah

s_M = [i for i in range(s_ah)]
s_N = [i for i in range(s_aw)]
s_sv = SV()
s_sv.set_A(s_A, [i for i in range(aw, s_aw)], s_N)
s_sv.set_x(s_x, s_N)
s_sv.set_B(np.identity(s_ah))
print(simplex(s_A, b, s_c, s_sv, s_M, s_N))
print(s_sv.x)


sv = SV()
sv.A = s_sv.A
sv.B = s_sv.B
sv.N = s_sv.N
sv.L = list_diff(N, sv.N)
sv.x = s_sv.x[N]
sv.N_plus = [i for i in N if sv.x[i] > 0]
sv.N_zero = [i for i in N if sv.x[i] == 0]
print(simplex(A, b, c, sv, M, N))


# check section
print(A.dot(sv.x))

print(sv.x)
# print(c)
# for i in N:
#     print(sv.x[i])
#     print(c[i])
    # print(sv.x[i] * c[i])
print(-sv.x.dot(c))











