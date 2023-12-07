import numpy as np
import matplotlib.pyplot as plt
import math

d = np.loadtxt("time_required.out", delimiter=";")

plt.rcParams.update({"font.size": 30})

# f = plt.figure()
f, ax = plt.subplots()

ax.plot(np.log(d[:, 0]), np.log(d[:, 1]), "k", linewidth=4.0)
ax.plot(np.log(d[:, 0]), np.log(d[:, 2]), "r", linewidth=4.0)
ax.plot(np.log(d[0:4, 0]), np.log(d[0:4, 3]), "b", linewidth=4.0)
ax.plot(np.log(d[0:4, 0]), np.log(d[0:4, 4]), "g", linewidth=4.0)

a, b = np.polyfit(np.log(d[1:6, 0]), np.log(d[1:6, 1]), 1)
print(a)
a1, b1 = np.polyfit(np.log(d[1:6, 0]), np.log(d[1:6, 2]), 1)
print(a1)
a2, b2 = np.polyfit(np.log(d[0:4, 0]), np.log(d[0:4, 3]), 1)
print(a2)
a3, b3 = np.polyfit(np.log(d[1:4, 0]), np.log(d[1:4, 4]), 1)
print(a3)
ax.plot(np.log(d[:, 0]), a * np.log(d[:, 0]) + b, "k--", linewidth=4.0)
ax.plot(np.log(d[:, 0]), a1 * np.log(d[:, 0]) + b1, "r--", linewidth=4.0)
ax.plot(np.log(d[0:4, 0]), a2 * np.log(d[0:4, 0]) + b2, "b--", linewidth=4.0)
ax.plot(np.log(d[0:4, 0]), a3 * np.log(d[0:4, 0]) + b3, "g--", linewidth=4.0)
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel("ln(frequency)")
ax.set_ylabel("ln(time)")
ax.legend(["Old IDSM", "New IDSM", "LUDSM"])
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off

# ax.set_xticks([1, 2, 3, 4], minor=True)
# ax.set_xticklabels([1, 2, 3, 4], minor=True)

plt.show()
