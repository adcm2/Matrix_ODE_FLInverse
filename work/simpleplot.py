import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("frawspec.out", delimiter=";")
dp = np.loadtxt("drawspec")
plt.plot(d[:, 0], d[:, 3], "k")
plt.plot(dp[:, 0], dp[:, 3], "r")

plt.xlabel("Frequency (mHz)")
plt.show()
print(len(dp))