import numpy as np
import matplotlib.pyplot as plt
import math

d = np.loadtxt("spectra_test.STAT1.Z")
dp = np.loadtxt("fspectra.r1.out.q4", delimiter=";")
dpp = np.loadtxt("fspectra.fullcouple.r1.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.N")
# dp = np.loadtxt("fspectra.r2.out.q4", delimiter=";")
# dpp = np.loadtxt("fspectra.fullcouple.r2.out.q4", delimiter=";")
# d = np.loadtxt("spectra_test.STAT1.E")
# dp = np.loadtxt("fspectra.r3.out.q4", delimiter=";")
# dpp = np.loadtxt("fspectra.fullcouple.r3.out.q4", delimiter=";")


plt.rcParams.update({"font.size": 30})

# f = plt.figure()
f, ax = plt.subplots(1,1)



mynum = max(abs(d[:,3]))
# # ax[0,0].subplot(2,1,1)
# ax[0].plot(d[:, 0], d[:, 3], "k", linewidth = 2.0)
# ax[0].plot(dp[:, 0], dp[:, 3], "r", linestyle = 'dashed', linewidth = 2.0)
# ax[0].plot(dpp[:, 0], 2*dpp[:, 3], "b", linestyle = 'dotted', linewidth = 2.0)
# ax[0].legend(['Old code', 'New code', "Full coupling"])

# ax[0].set_ylabel('Arb. units')
# ax[0].set_yticks([])
# ax[0].set_ylim(0,mynum)

# # # relative difference
dout = np.zeros(len(d))
dout2 = np.zeros(len(d))

for i in range(0,len(d)-1,1):
    dout[i] = abs(2*dpp[i, 3] - d[i,3])/mynum
    dout2[i] = abs(2*dpp[i, 3] - dp[i,3])/mynum

# # plt.subplot(2,1,2)
# ax[1].plot(d[:, 0], dout, "b", linewidth = 2.0)
# ax[1].plot(d[:, 0], dout2, "r", linewidth = 2.0)
# ax[1].ticklabel_format(axis='both', style='sci', scilimits=(4,4))
# ax[1].set_xlabel('Frequency (mHz)')
# ax[1].set_ylabel('% difference (x$ 10^{-5}$)')
# ax[1].legend(['Difference between full coupling and old', 'Difference between full coupling and new'])
# # t = ax[1].yaxis.get_offset_text()
# ax[1].yaxis.get_offset_text().set_visible(False)
# # t.set_x(-0.1)
# # text(x, y, "1e4")
# # plt.ticklabel_format(axis = 'y', style = 'sci')

ax.plot(d[:, 0], dout, "b", linewidth = 2.0)
ax.plot(d[:, 0], dout2, "r", linewidth = 2.0)
ax.ticklabel_format(axis='both', style='sci', scilimits=(4,4))
ax.set_xlabel('Frequency (mHz)')
ax.set_ylabel('% difference (x$ 10^{-5}$)')
ax.legend(['Difference between full coupling and old', 'Difference between full coupling and new'])
# t = ax[1].yaxis.get_offset_text()
# ax.yaxis.get_offset_text().set_visible(False)
# t.set_x(-0.1)
# text(x, y, "1e4")
# plt.ticklabel_format(axis = 'y', style = 'sci')

plt.show()
