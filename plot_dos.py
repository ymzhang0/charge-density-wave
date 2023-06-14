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


def plot_dos(DOS_data, FERMI_LEVEL):


    plt.rc('font', family='Times New Roman')


    E = DOS_data[:, 0]   - FERMI_LEVEL
    DE = DOS_data[:, 1]


    fig, ax = plt.subplots(figsize=(6, 12))


    ax.set_ylim(-2, 2)

    ax = plt.gca();  # 获得坐标轴的句柄
    ax.spines['bottom'].set_linewidth(2);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(2);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(2);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(2);  ####设置上部坐标轴的粗细
    ax.tick_params(axis="x", direction="in")
    ax.tick_params(axis="y", direction="in")




    # plt.xticks([-2, -1, 0, 1, 2])
    plt.xticks([])
    plt.yticks([])
    plt.axline((0, 0), (np.max(DE), 0), linewidth=1.5, color='grey')
    # ax.set_yticklabels(["-4", "-2", "0", "2", "4"],size=32,position=(-0.02,0), fontdict={'family': 'Arial','weight': 'bold'})
    # ax.set_xticklabels(["-2", "-1", "0", "1", "2"], size=20, position=(-0.02, 0),
    #                    fontdict={'family': 'Arial', 'weight': 'bold'})


    a = plt.plot(DE, E, color='red', alpha=0.8, lw=2)

    # ax.set_xlabel("Energy(eV)", size=24, fontdict={'family': 'Arial', 'weight': 'bold'})
    # ax.set_ylabel("Density of States", size=24, fontdict={'family': 'Arial', 'weight': 'bold'})
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.ioff()
    plt.savefig(root + "//DOS/pdf" + "//%s%s2.pdf" %(cation, anion), format='pdf', dpi=2000)



ifile = 0
print("Presen file: %s" % root)
for cation in cations:
    for anion in anions:
        EMIN, EMAX, NSTEPS, FERMI_ENERGY = 0.0, 0.0, 0, 0.0
        DOS = []
        with open(root + "\\DOS\\%02d_%s%s2\\DOSCAR" %(ifile, cation, anion), "r") as f:
            lines = f.readlines()
            EMIN, EMAX, NSTEPS, FERMI_ENERGY = lines[5].strip().split()[:4]
            for line in lines[7: 6+int(NSTEPS)]:
                DOS.append(line.strip().split())
        DOS = np.array(DOS).astype(float)

        plot_dos(DOS, float(FERMI_ENERGY))




        ifile += 1
