import numpy as np
from math import sqrt, fabs, sin, cos, log
import matplotlib.pyplot as plot
import matplotlib.patches as mpatches


def estimation_gr(E):
    return log(E, 2 / (1 + sqrt(5)))


def estimation_us(E):
    return 3 * log(E, 1/3)


E = 1
legends = []
for i in range(8):
    plot.plot(i, estimation_gr(E), 'r*')
    plot.plot(i, estimation_us(E), 'g*')
    E /= 10
plot.xlabel('i')
plot.ylabel('calls')
legends.append(mpatches.Patch(color='r', label='золотое сечение'))
legends.append(mpatches.Patch(color='g', label='равномерный поиск'))
plot.legend(handles=legends)
plot.show()