import numpy as np
import matplotlib.pyplot as plt
import math


d = np.loadtxt("time_required.out", delimiter=";")

plt.rcParams.update({"font.size": 24})

f = plt.figure()
plt.yscale("log")
plt.plot(d[:, 0], d[:, 1], "k")
plt.plot(d[:, 0], d[:, 2], "r")
plt.legend(['Old code', 'New code'])
plt.show()