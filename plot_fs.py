import os
from math import sqrt, pi
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import rcParams
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.patches as mpatches

cations = ["Ti", "V", "Nb", "Ta"]
anions = ["S", "Se", "Te"]
phases = ["2H", "1T"]
SYMMETRY_ops = [[[1, 0, 0],
                 [0, 1, 0],
                 [0, 0, 1]
                 ],
                [[0, 1, 0],
                 [1, 0, 0],
                 [0, 0, 1]
                 ],
                [[1, 1, 0],
                 [0, -1, 0],
                 [0, 0, 1]
                 ],
                [[0, -1, 0],
                 [1, 1, 0],
                 [0, 0, 1]],
                [[1, 1, 0],
                 [-1, 0, 0],
                 [0, 0, 1]],
                [[1, 0, 0],
                 [-1, -1, 0],
                 [0, 0, 1]],
                [[-1, 0, 0],
                 [0, -1, 0],
                 [0, 0, 1]
                 ],
                [[0, -1, 0],
                 [-1, 0, 0],
                 [0, 0, 1]
                 ],
                [[-1, -1, 0],
                 [0, 1, 0],
                 [0, 0, 1]
                 ],
                [[0, 1, 0],
                 [-1, -1, 0],
                 [0, 0, 1]],
                [[-1, -1, 0],
                 [1, 0, 0],
                 [0, 0, 1]],
                [[-1, 0, 0],
                 [1, 1, 0],
                 [0, 0, 1]]
                ]
root = os.getcwd()

fs_dat_dir = root + "\\fs_ibz"
fs_fig_dir = root + "\\fs_xy"
RECIPROCAL_LATTICE = np.array([[sqrt(3)/2  ,  1/2        ,  0.000000000],
                               [0.000000000,  1.000000000,  0.000000000],
                               [0.000000000, -0.000000000,  1.000000000]
                              ])

def load_data(localdir):
    KPOINTS_frac = []
    ENERGIES = []
    OCCUPATIONS = []
    NELEC, NKPOINTS, NBANDS, FERMI_ENERGY = 0, 0, 0, 0
    with open(root + "\\fs_ibz_kzhalf" +localdir + "DOSCAR", 'r') as f:
        line = f.readlines()[5].strip().split()
        FERMI_ENERGY = float(line[-2])

    with open(root + "\\fs_ibz_kzhalf" +localdir + "EIGENVAL", 'r') as f:
        line = f.readlines()[5].strip().split()
        NELEC, NKPOINTS, NBANDS = [int(i) for i in line]

    with open(root + "\\fs_ibz_kzhalf" +localdir + "EIGENVAL", 'r') as f:
        lines = f.readlines()[7:]
        for line in lines:
            line = line.strip().split()
            if len(line) == 4:
                KPOINTS_frac.append(line)
            if len(line) == 3:
                ENERGIES.append(line[1])
                OCCUPATIONS.append(line[2])
    print("fermi energy for %02d_%s%s2  is %.4f" %(ifile, cation, anion, FERMI_ENERGY))
    return NELEC, NKPOINTS, NBANDS, FERMI_ENERGY, KPOINTS_frac, ENERGIES, OCCUPATIONS

ifile = 0
print("Presen file: %s" % root)

for cation in cations:
    for anion in anions:
# for cation in ["V"]:
#     for anion in ["S"]:
        FERMI_BAND = []
        NELEC, NKPOINTS, NBANDS, FERMI_ENERGY, KPOINTS_frac, ENERGIES, OCCUPATIONS = load_data("\\%02d_%s%s2\\" % (ifile, cation, anion))

        KPOINTS_frac = np.array(KPOINTS_frac).reshape(NKPOINTS, 4).astype(float)
        ENERGIES = np.array(ENERGIES).reshape(NKPOINTS, NBANDS).astype(float) - FERMI_ENERGY
        OCCUPATIONS = np.array(OCCUPATIONS).reshape(NKPOINTS, NBANDS).astype(float)


        KPOINTS_frac_fbz = np.array([])


        for isym in SYMMETRY_ops:
            for ik in KPOINTS_frac:
                k_eq = np.dot(np.array(isym), ik[:3])
                KPOINTS_frac_fbz = np.append(KPOINTS_frac_fbz, k_eq)

        KPOINTS_frac_fbz = KPOINTS_frac_fbz.reshape(NKPOINTS * 12, 3)
        KPOINTS_cart_fbz = np.dot(KPOINTS_frac_fbz[:, :3], RECIPROCAL_LATTICE)

        for iband in range(NBANDS):
            lower_bound, upper_bound = np.min(ENERGIES[:, iband]), np.max(ENERGIES[:, iband])
            if (lower_bound < 0) & (upper_bound > 0):
                FERMI_BAND.append([iband, lower_bound, upper_bound])


        fig, ax = plt.subplots(figsize=(6, 6))
        frame = plt.gca()

        frame.axes.get_xaxis().set_visible(False)
        frame.axes.get_yaxis().set_visible(False)

        for i in ['top', 'right', 'bottom', 'left']:
            ax.spines[i].set_visible(False)

        plt.tick_params(bottom=False, top=False, left=False, right=False)

        ax.set_ylim(-0.8, 0.8)
        ax.set_xlim(-0.8, 0.8)
        allbands = []
        for ifermiband in FERMI_BAND:
            allbands.append(ifermiband[0])
            dis_energy = np.exp(-ENERGIES[:, ifermiband[0]] ** 2 / 0.001)
            a = plt.scatter(KPOINTS_cart_fbz[:, 0], KPOINTS_cart_fbz[:, 1], s=2, c=np.tile(dis_energy, (12, 1)),
                            cmap="Blues")
            # a=plt.scatter(KPOINTS_cart_fbz[:, 0],KPOINTS_cart_fbz[:, 1], s=1)
            plt.ioff()
            plt.savefig(root + "\\fs_ibz_kzhalf\\pdf" + "\\%s%s2_band%d.pdf" % (cation, anion, ifermiband[0]), format='pdf', dpi=2000)

        G = [0, 0, 0]
        M = np.dot(np.array([1/2, 0, 0]), RECIPROCAL_LATTICE)
        K = np.dot(np.array([1 /3, 1 /3, 0]), RECIPROCAL_LATTICE)


        dis_energy = np.sum(np.exp(-ENERGIES[:, allbands] ** 2 / 0.005), axis=1)
        a = plt.scatter(KPOINTS_cart_fbz[:, 0], KPOINTS_cart_fbz[:, 1], s=2, c=np.tile(dis_energy, (12, 1)),
                        cmap="Blues")
        # a=plt.scatter(KPOINTS_cart_fbz[:, 0],KPOINTS_cart_fbz[:, 1], s=1)
        plt.ioff()
        plt.savefig(root + "\\fs_ibz_kzhalf\\pdf" + "\\%s%s2_sum.pdf" % (cation, anion), format='pdf', dpi=2000)

        # plot_fs(cation, anion)

        ifile += 1
