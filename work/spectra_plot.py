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

plt.rcParams.update({"font.size": 24})

f = plt.figure()

# print(d[0:3,:])
# print(dp[0:3,:])
mynum = max(abs(d[:,3]))
plt.subplot(2,1,1)
plt.plot(d[:, 0], d[:, 3], "k")
plt.plot(dp[:, 0], dp[:, 3], "r")
plt.plot(dpp[:, 0], dpp[:, 3], "b")
plt.legend(['Fortran', 'C++'])

plt.ylabel('Comparison')
# plt.ylim(0,mynum)


# # # relative difference
dout = np.zeros(len(d))

for i in range(0,len(d)-1,1):
    dout[i] = abs(dp[i, 3] - d[i,3])/mynum

plt.subplot(2,1,2)
plt.plot(d[:, 0], dout, "k")
plt.xlabel('Frequency (mHz)')
plt.ylabel('% difference')



plt.show()
