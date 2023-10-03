import numpy as np
import matplotlib.pyplot as plt
import math

d = np.loadtxt("spectra_test.STAT1.Z")
dp = np.loadtxt("fspectra.r1.out", delimiter=";")

plt.rcParams.update({"font.size": 24})

f = plt.figure()

# plt.subplot(3,1,1)
# plt.plot(d[:, 0], d[:, 3], "k")
# plt.plot(dp[:, 0], dp[:, 3], "r")
# plt.legend(['Fortran', 'C++'])
# plt.title('Direct comparison of two solvers')

# plt.subplot(3,1,2)
# plt.plot(d[:, 0], abs(dp[:, 3] - d[:,3]), "r")
# plt.title('Absolute difference')

# # relative difference
# dout = np.zeros(len(d))
# mynum = max(abs(d[:,3]))
# for i in range(0,len(d)-1,1):
#     # dout[i-1] = abs(dp[i-1, 3] - d[i,3])/d[i,3]*100*(1 - math.exp(-(100 * abs(d[i,3])/mynum)**2))
#     if abs(d[i,3]) > 0.001 * mynum:
#         dout[i] = abs(dp[i, 3] - d[i,3])/d[i,3]*100

# plt.subplot(3,1,3)
# plt.plot(d[:, 0], dout, "r")
# plt.title('Relative difference')

plt.subplot(2,1,1)
plt.plot(d[:, 0], d[:, 3], "k")
plt.plot(dp[:, 0], dp[:, 3], "r")
plt.legend(['Fortran', 'C++'])
plt.title('Direct comparison of two solvers')


# relative difference
dout = np.zeros(len(d))
mynum = max(abs(d[:,3]))
for i in range(0,len(d)-1,1):
    # dout[i-1] = abs(dp[i-1, 3] - d[i,3])/d[i,3]*100*(1 - math.exp(-(100 * abs(d[i,3])/mynum)**2))
    if abs(d[i,3]) > 0.01 * mynum:
        dout[i] = abs(dp[i, 3] - d[i,3])/d[i,3]*100

plt.subplot(2,1,2)
plt.plot(d[:, 0], dout, "r")
plt.title('Percentage difference')

# # # plt.plot(d[1:len(d), 0], abs(dp[:, 3] - d[1:len(d),3])/d[1:len(d),3]*100*(1 - math.exp(-abs(d[1:len(d),3]))), "r")
# plt.plot(d[1:len(d), 0], dout, "r")
# # plt.title('"Relative" difference')



# relative difference
# plt.subplot(3,1,3)
# plt.plot(d[:, 0], abs((dp[:, 3] - d[:,3])), "r")
# plt.subplot(2,1,1)

# plt.plot(d[:, 0], d[:, 3], "k")
# plt.plot(dp[:, 0], dp[:, 3], "r")
# plt.legend(['Fortran', 'C++'])
# plt.title('Direct comparison of two solvers')

# dout = np.zeros(len(d)-1)
# # print(len(dout))
# plt.subplot(3,1,2)
# mynum = max(abs(d[:,3]))
# for i in range(1,len(d),1):
#     # dout[i-1] = abs(dp[i-1, 3] - d[i,3])/d[i,3]*100*(1 - math.exp(-(100 * abs(d[i,3])/mynum)**2))
#     if abs(d[i,3]) > 0.001 * mynum:
#         dout[i-1] = abs(dp[i-1, 3] - d[i,3])/d[i,3]*100

# # # plt.plot(d[1:len(d), 0], abs(dp[:, 3] - d[1:len(d),3])/d[1:len(d),3]*100*(1 - math.exp(-abs(d[1:len(d),3]))), "r")
# plt.plot(d[1:len(d), 0], dout, "r")
# plt.title('"Relative" difference')



# plt.subplot(2,1,2)
# plt.plot(d[1:len(d), 0], abs(dp[:, 3] - d[1:len(d),3]), "r")




# plt.plot(d[1:len(d), 0], dout, "r")
# plt.title('Absolute difference')

# print(len(d[1:len(d),3]))
# print(len(dout))
# print(d[1:10,0])
# print(dp[0:10,0])

# print(d[9435:9442,0])
# print(dp[9435:9441,0])

# plt.gca().set_yticks([])
# plt.yticks([])
# plt.xlabel("Frequency (mHz)")
plt.show()
