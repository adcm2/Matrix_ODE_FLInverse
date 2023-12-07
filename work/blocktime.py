import numpy as np
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("spectra_test.STAT1.Z")
d = np.loadtxt("time_block.out", delimiter=";")
# dpp = np.loadtxt("fspectra.fullcouple.r1.out.q4", delimiter=";")
# dp = np.loadtxt("fspectra.r1.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.N")
# dp = np.loadtxt("fspectra.r2.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.E")
# dp = np.loadtxt("fspectra.r3.out.q4", delimiter=";")

plt.rcParams.update({"font.size": 30})

f = plt.figure()


plt.plot(d[:, 0], d[:, 1], "k", linewidth=4.0)
plt.xlabel("Width of target block (mHz)")

plt.ylabel("Time (s)")
# plt.ylim(0,mynum)


# # # relative difference
# dout = np.zeros(len(d))
# dout2 = np.zeros(len(d))

# for i in range(0, len(d) - 1, 1):
#     dout[i] = abs(2*dpp[i, 3] - d[i, 3]) / mynum
#     dout2[i] = abs(2*dpp[i, 3] - dp[i, 3]) / mynum

# plt.subplot(2, 1, 2)
# plt.plot(d[:, 0], dout, "k")
# plt.plot(d[:, 0], dout2, "r")
# plt.xlabel("Frequency (mHz)")
# plt.ylabel("% difference")


plt.show()
