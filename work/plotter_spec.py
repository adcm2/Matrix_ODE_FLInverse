import numpy as np
import matplotlib.pyplot as plt

d = np.loadtxt("test.out", delimiter=";")
d2 = np.loadtxt("testrawspec.out")
d3 = np.loadtxt("fspectra.r1.out", delimiter=";")
d4 = np.loadtxt("spectra_test.STAT1.Z")

# print(len(d))
# print(len(d4))
idx = 0
idx2 = 0
# while idx < d2.shape[0]:
#     idx2 = 0
#     while idx2 < d2.shape[1]:
#         if abs(d2[idx,idx2]) > 1:
#             d2[idx,idx2] = 0
#         idx2 += 1
#     idx+=1
# plt.plot(1000*d[:,0],d[:,1],'k')
plt.plot(d3[:, 0], d3[:, 3], "k")
plt.plot(d4[:, 0], d4[:, 3] , "r")
# plt.plot(1000*d[:,0], d2[0:len(d),2],'r')
# plt.plot(1000*d[:,0],d[:,2],'r')
# plt.plot(1000*d[:,0],d[:,3],'b')

# # plt.plot(dppp[:,0],dppp[:,3],'r--')
# # plt.plot(dpppp[:,0],dpppp[:,3],'b--')

# plt.plot(d[:,0],d[:,3] - dppp[:,3],'k')
# plt.plot(d[:,0],d[:,3] - dpppp[:,3],'r')


# plt.gca().set_yticks([])
# plt.yticks([])
# plt.xlim(0,1)
plt.xlabel("Frequency (mHz)")
plt.show()
