import numpy as np
import matplotlib.pyplot as plt

all_energies_a = np.loadtxt("all_energies_in_outc")
#
# def derivitave(x, coeff):
#     y = np.array([])
#     der = np.array([])
#     for j in range(1, len(coeff)+1):
#         der = np.append(der, (j - 1) * coeff[j])
#
#     for i in range(len(x)):
#         temp = 0.0
#         for j in range(len(coeff)):
#             temp += coeff[j] * x[i] ** j
#         y = np.append(y, temp)
#
#     return
#


def derivitave(coeff):
    degree = len(coeff)
    der = np.array([])
    for j in range(0, degree):
        der = np.append(der, (degree - j - 1) * coeff[j])
    return der

relaxed_a = []

for imat in range(len(all_energies_a)):
    print("Finding energies minumum for mat %02d " % (imat))
    degree = 5
    x = np.linspace(all_energies_a[imat, 0], all_energies_a[imat, 1], 20)
    y = all_energies_a[imat, 2:]
    coeff = np.polyfit(x, y, degree)
    der_coeff = derivitave(coeff)
    p = np.poly1d(coeff)
    d = np.poly1d(der_coeff)
    plt.plot(x, d(x))
    for iroot in d.roots:
        if abs(iroot.imag) < 1e-10 and np.abs(iroot) < all_energies_a[imat, 1] and np.abs(iroot) > all_energies_a[imat, 0]:
            print(iroot.imag )
            print("Found energies minumum for mat %02d : %.6f + %.6fj" %(imat, iroot.real, iroot.imag))
            relaxed_a.append(iroot.real)

if len(relaxed_a) == 12:
    print("Done")
    
np.savetxt("Relaxed_c", np.array(relaxed_a))