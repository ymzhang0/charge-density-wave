import os
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import rcParams
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.patches as mpatches


cations = ["Ti", "V", "Nb", "Ta"]
anions = ["S", "Se", "Te"]
phases = ["2H", "1T"]

root = os.getcwd()

band_dat_dir = "\\xml_for_phonon"

ifile = 0
print("Presen file: %s" % root)

KTICKS = []
KLABELS = ["Γ", "M", "K", "Γ", "A", "L", "H", "A"]
mat = "05_VTe2"
dat = "VTe2_442.dat"
natom = 9
with open(".\\xml_for_phonon\\%s\\%s" %(mat, dat), "r") as f:
    lines = f.readlines()
    for line in lines[1].strip().split()[2:]:
        if line == []:
            break
        else:
            KTICKS.append(float(line))


print(KTICKS)
plt.rc('font', family='Times New Roman')

band = np.loadtxt(".\\xml_for_phonon\\%s\\%s" %(mat, dat))

k = band[:, 0].reshape(natom, len(band)//natom)
E = band[:, 1].reshape(natom, len(band)//natom)


fig, ax = plt.subplots(figsize=(6, 6))


ax = plt.gca();
ax.spines['bottom'].set_linewidth(2);  ###设置底部坐标轴的粗细
ax.spines['left'].set_linewidth(2);  ####设置左边坐标轴的粗细
ax.spines['right'].set_linewidth(2);  ###设置右边坐标轴的粗细
ax.spines['top'].set_linewidth(2);  ####设置上部坐标轴的粗细
ax.tick_params(axis="x", direction="in")
ax.tick_params(axis="y", direction="in")



ax.set_xticks(KTICKS)
ax.set_xticklabels(KLABELS, size=20, position=(0, -0.02),
                   fontdict={'family': 'Arial', 'weight': 'bold'})

# plt.yticks([-2, -1, 0, 1, 2])
# # ax.set_yticklabels(["-4", "-2", "0", "2", "4"],size=32,position=(-0.02,0), fontdict={'family': 'Arial','weight': 'bold'})
ax.set_yticklabels([], size=20, position=(-0.02, 0),
                   fontdict={'family': 'Arial', 'weight': 'bold'})
for xtick in KTICKS[1:-1]:
    ax.axvline(x=xtick, ls="--", c="black", alpha=0.6, linewidth=2)

# ax.set_ylim(-2, 2)
ax.set_xlim(0, KTICKS[-1])
# a = plt.plot(k, E, color='black', alpha=0.8, lw=1)
plt.axline((0, 0), (KTICKS[-1], 0), linewidth=1, color='Grey')
for i in range(natom):
    plt.plot(k[i], E[i], color='black', alpha=0.8, lw=1)

a = plt.scatter(k, E, s=3, facecolors='none', edgecolors='red', alpha=0.8)
ax.set_ylabel("Phonon Frequency", size=20, fontdict={'family': 'Arial', 'weight': 'bold'})
plt.ioff()
plt.savefig(".\\xml_for_phonon\\pdf\\%s.pdf" %(mat) , format='pdf', dpi=2000)

ifile += 1
