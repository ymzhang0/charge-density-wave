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




ifile = 0
print("Presen file: %s" % root)
for cation in cations:
    for anion in anions:
        os.chdir(root + "\\band\\%02d_%s%s2" %(ifile, cation, anion))

        KLABELS = []

        with open("KLABELS", "r") as f:

            lines = f.readlines()
            for line in lines[1:]:
                line = line.strip().split()
                if line != []:
                    KLABELS.append(line)
                else:
                    break

        print(KLABELS)
        plt.rc('font', family='Times New Roman')

        FERMI_ENERGY = -2.473595
        cation_band = np.loadtxt("PBAND_%s.dat" % cation)
        anion_band = np.loadtxt("PBAND_%s.dat" % anion)

        k = cation_band[:, 0]
        E = cation_band[:, 1]  # - FERMI_ENERGY

        cation_s = cation_band[:, 2]
        cation_p = cation_band[:, 3:6]
        cation_d = cation_band[:, 6:11]

        anion_s = anion_band[:, 2]
        anion_p = anion_band[:, 3:6]
        anion_d = anion_band[:, 6:11]

        fig, ax = plt.subplots(figsize=(6, 6))

        # ax.set_ylim(-2, 2)
        # ax.set_xlim(0, 2.785)

        ax = plt.gca();  # 获得坐标轴的句柄
        ax.spines['bottom'].set_linewidth(2);  ###设置底部坐标轴的粗细
        ax.spines['left'].set_linewidth(2);  ####设置左边坐标轴的粗细
        ax.spines['right'].set_linewidth(2);  ###设置右边坐标轴的粗细
        ax.spines['top'].set_linewidth(2);  ####设置上部坐标轴的粗细
        ax.tick_params(axis="x", direction="in")
        ax.tick_params(axis="y", direction="in")

        # ax.arrow(2.5262, -0.6, 0, 1.92,
        #              width=0.01,
        #              length_includes_head=True,
        #               head_width=0.05,
        #               head_length=0.1,
        #              fc='r',
        #              ec='r')

        xticks = [float(i[1]) for i in KLABELS]
        xlabels = ["Γ" if i[0] == "GAMMA" else i[0] for i in KLABELS]
        xlabels = ["A" if i[0] == "A|L" else i[0] for i in xlabels]

        ax.set_xticks(xticks[:-2])
        ax.set_xticklabels(xlabels[:-2], size=20, position=(0, -0.02),
                           fontdict={'family': 'Arial', 'weight': 'bold'})

        plt.yticks([-2, -1, 0, 1, 2])
        # ax.set_yticklabels(["-4", "-2", "0", "2", "4"],size=32,position=(-0.02,0), fontdict={'family': 'Arial','weight': 'bold'})
        ax.set_yticklabels(["-2", "-1", "0", "1", "2"], size=20, position=(-0.02, 0),
                           fontdict={'family': 'Arial', 'weight': 'bold'})
        for xtick in xticks[:-3]:
            ax.axvline(x=xtick, ls="--", c="black", alpha=0.6, linewidth=2)

        # plt.axline((0, 0), (2.785, 0), linewidth=2, color='black')

        # ax.axvline(x=1.80156577,ls="-",c="black",alpha=0.6, linewidth=2)
        # ax.axvline(x=1.75810,ls="-",c="black",alpha=0.6, linewidth=8)

        # a = plt.plot(x, y, c='black', alpha=1)
        # b = plt.plot(x, y, c='r', alpha=1)
        ax.set_ylim(-2, 2)
        ax.set_xlim(0, xticks[-3])
        a = plt.plot(k, E, color='black', alpha=0.8, lw=1)

        a = plt.scatter(k, E, s=np.sum(cation_d, axis=1) * 20, facecolors='none', edgecolors='r', alpha=0.8)
        a = plt.scatter(k, E, s=np.sum(anion_p, axis=1) * 20, facecolors='none', edgecolors='green', alpha=0.8)
        ax.set_ylabel("Energy(eV)", size=24, fontdict={'family': 'Arial', 'weight': 'bold'})
        plt.ioff()
        plt.savefig(root + "\\band\\pdf\\%02d_%s%s2.pdf" % (ifile, cation, anion), format='pdf', dpi=2000)

        ifile += 1
