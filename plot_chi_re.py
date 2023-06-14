import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt

from matplotlib.pyplot import rcParams
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.patches as mpatches
import yaml
import os

cations = ["Ti", "V", "Nb", "Ta"]
anions = ["S", "Se", "Te"]
phases = ["2H", "1T"]

# plt.rc('font',family='Times New Roman')
Hartree2eV = 27.211386245988
Boltzmann_eV = 8.617333262145 * 1e-5 # eV
Bohr_radius = 0.52917721067 # ang

Nkx = Nky = 105
Nqx = Nqy = 105

current_dir = os.getcwd()
chi_dir = current_dir + "\\chi_re_ip"
to_fig_dir = current_dir + "\\chi_re_plot"

RECIPROCAL_LATTICE = np.array([[sqrt(3)/2  ,  1/2        ,  0.000000000],
                               [0.000000000,  1.000000000,  0.000000000],
                               [0.000000000, -0.000000000,  1.000000000]
                              ])
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
# x, y, z = np.meshgrid(np.linspace(0, (Nqx - 1) / Nqx, Nqx), np.linspace(0, (Nqy - 1) / Nqy, Nqy), [0])
# QPOINTS_frac = np.array([x.ravel(), y.ravel(), z.ravel()]).T
# QPOINTS_cart = np.dot(QPOINTS_frac, RECIPROCAL_LATTICE)
#
# QPOINTS_cart_extended = np.vstack([QPOINTS_cart, QPOINTS_cart - RECIPROCAL_LATTICE[0, :],
#                                        QPOINTS_cart - RECIPROCAL_LATTICE[1, :], QPOINTS_cart - RECIPROCAL_LATTICE[0, :] - RECIPROCAL_LATTICE[1, :]])

QPOINTS_frac = np.loadtxt(".\\chi_ip\\KPOINTS_63_63_1", skiprows=3)

ifile = 0
for cation in cations:
    for anion in anions:
        QPOINTS_frac_fbz = np.array([])
        chi_q = np.loadtxt(".\\chi_ip\\%02d_%s%s2" %(ifile, cation, anion), skiprows=1)
        chi_q_re = np.sum(chi_q[:, 0:chi_q.shape[1]:2], axis=1)
        chi_q_im = np.sum(chi_q[:, 1:chi_q.shape[1]:2], axis=1)

        # for iband in range(len(fermiband)):
        #     chi_extended = np.tile(chi_q[iband].flatten(), (4, 1))
        #     fig,ax=plt.subplots(figsize=(6,6))
        #
        #     for i in ['top', 'right', 'bottom', 'left']:
        #         ax.spines[i].set_visible(False)
        #
        #     ax.set_ylim(-1.5, 1.5)
        #     ax.set_xlim(-1.5, 1.5)
        #     a=plt.scatter(QPOINTS_cart_extended[:, 0],QPOINTS_cart_extended[:, 1], s=1, c=chi_extended,  cmap="Blues")
        #     plt.ioff()
        #     plt.savefig(to_fig_dir + "\\%02d_%s%s2_band%d.pdf" %(ifile, cation, anion, int(fermiband[iband])),format='pdf', dpi=2000)
        for isym in SYMMETRY_ops:
            for ik in QPOINTS_frac:
                k_eq = np.dot(np.array(isym), ik[:3])
                QPOINTS_frac_fbz = np.append(QPOINTS_frac_fbz, k_eq)
        QPOINTS_frac_fbz = QPOINTS_frac_fbz.reshape(QPOINTS_frac.shape[0] * 12, 3)
        QPOINTS_cart_fbz = np.dot(QPOINTS_frac_fbz[:, :3], RECIPROCAL_LATTICE)
        chi_re_fbz = np.tile(chi_q_re, (12, 1))
        chi_im_fbz = np.tile(chi_q_im, (12, 1))

        fig,ax=plt.subplots(figsize=(6,6))

        for i in ['top', 'right', 'bottom', 'left']:
            ax.spines[i].set_visible(False)

        ax.set_ylim(-1.5, 1.5)
        ax.set_xlim(-1.5, 1.5)
        a=plt.scatter(QPOINTS_cart_fbz[:, 0],QPOINTS_cart_fbz[:, 1], s=3, c=-chi_re_fbz,  cmap="Blues")
        plt.ioff()
        plt.savefig(".\\chi_ip\\pdf\\%s%s2_re.pdf" %(cation, anion),format='pdf', dpi=2000)
        print("%s%s2 done!" %(cation, anion))
        ifile += 1