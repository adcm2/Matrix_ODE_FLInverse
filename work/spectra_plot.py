import numpy as np
import matplotlib.pyplot as plt
import math

d = np.loadtxt("spectra_test.STAT1.Z")
dp = np.loadtxt("fspectra.r1.out.q4", delimiter=";")
dpp = np.loadtxt("fspectra.fullcouple.r1.out.q4", delimiter=";")
# dp = np.loadtxt("fspectra.r1.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.N")
# dp = np.loadtxt("fspectra.r2.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.E")
# dp = np.loadtxt("fspectra.r3.out.q4", delimiter=";")

plt.rcParams.update({"font.size": 30})

# f = plt.figure()
f,ax = plt.subplots()

# print(d[0:3,0])
# print(dp[0:3,0])
# print(dpp[0:3,0])
mynum = max(abs(d[:, 3]))
# plt.subplot(2, 1, 1)
# plt.plot(d[:, 0], d[:, 3], "k")
# plt.plot(dp[:, 0], dp[:, 3], "r")
# plt.plot(dpp[:, 0], dpp[:, 3] * 2, "b")
# plt.legend(["Fortran", "C++"])

# plt.ylabel("Comparison")
# plt.ylim(0,mynum)


# # # relative difference
dout = np.zeros(len(d))
dout2 = np.zeros(len(d))

for i in range(0, len(d) - 1, 1):
    dout[i] = abs(2*dpp[i, 3] - d[i, 3]) / mynum
    dout2[i] = abs(2*dpp[i, 3] - dp[i, 3]) / mynum

# plt.subplot(2, 1, 2)
ax.plot(d[:, 0], dout, "k", linewidth = 2.0)
ax.plot(d[:, 0], dout2, "r", linewidth = 2.0)
ax.set_xlabel("Frequency (mHz)")
ax.set_ylabel("% difference (x$10^{-5}$)")
ax.legend(["Difference between eigensolution and old IDSM code", "Difference between eigensolution and new IDSM code"])
ax.ticklabel_format(axis = 'y', style = 'sci', scilimits = (0,0))

plt.show()
