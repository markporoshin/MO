import numpy as np


class SV:
    def set_A(self, A, N, M):
        self.A = A[:,N]
        self.N = N.copy()
        self.L = list_diff(M, N)


    def set_B(self, B):
        self.B = B.copy()


    def set_x(self, x, N):
        self.x = x
        self.N_plus = [i for i in N if x[i] > 0]
        self.N_zero = [i for i in N if x[i] == 0]


def list_diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif