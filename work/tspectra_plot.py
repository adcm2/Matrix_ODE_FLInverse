import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("tspectra.r1.out", delimiter=";")
d2 = np.loadtxt('time_test.STAT1.Z')

# Initialise the subplot function using number of rows and columns
# figure, axis = plt.subplots(1, 2)
plt.plot(d[:, 0], d[:, 1], "k")
plt.plot(d2[:,0],d2[:,1],'r')

# plt.ylim(-0.00001, 0.00001)
plt.show()
