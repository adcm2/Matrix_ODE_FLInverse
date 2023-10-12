import numpy as np
import matplotlib.pyplot as plt
import math

# d = np.loadtxt("spectra_test.STAT1.Z")
# dp = np.loadtxt("fspectra.r1.out", delimiter=";")
d = np.loadtxt("fspectra.r1.out.q1", delimiter=";")
dp = np.loadtxt("fspectra.r1.out.q4", delimiter=";")
dpp = np.loadtxt("fspectra.r1.out.q8", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.N")
# dp = np.loadtxt("fspectra.r2.out", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.E")
# dp = np.loadtxt("fspectra.r3.out", delimiter=";")

plt.rcParams.update({"font.size": 24})

f = plt.figure()
# print(len(d))
# print(len(dp))
# print(len(dpp))
# print(d[0:10,0])
# print(dp[0:7,0])
# print(dpp[0:7,0])
# print(d[len(d)-4:len(d)-1,0])
# print(dp[len(dp)-1])
# print(dpp[len(dpp) -1])

mynum = max(abs(d[:,3]))
plt.subplot(2,1,1)
plt.plot(d[:, 0], d[:, 3], "k")
plt.plot(dp[:, 0], dp[:, 3], "r")
plt.plot(dpp[:, 0], dpp[:, 3], "b")
plt.legend(['q = 1', 'q = 4', 'q = 8'])

plt.ylabel('Comparison')
# plt.ylim(0,mynum)

# relative difference
dout = np.zeros(len(dpp))
dout2 = np.zeros(len(dpp))

for i in range(0,len(dpp)-1,1):
    dout[i] = abs(dpp[i, 3] - d[i+7,3])/mynum

for i in range(0,len(dpp)-1,1):
    dout2[i] = abs(dp[i+1,3] - dpp[i,3])/mynum

plt.subplot(2,1,2)
plt.plot(dpp[:, 0], dout, "k")
plt.plot(dpp[:, 0], dout2, "b")
plt.legend(['q = 1', 'q = 4'])
plt.xlabel('Frequency (mHz)')
plt.ylabel('% error')



plt.show()
