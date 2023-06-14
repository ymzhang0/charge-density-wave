import os
from math import sqrt

cations = ["Ti", "V", "Nb", "Ta"]
anions = ["S", "Se", "Te"]
phases = ["2H", "1T"]

root = os.getcwd()

reference_alat_tmd = {"2D_2H_TiS2": [3.348], "2D_1T_TiS2": [3.417], "bulk_2H_TiS2": [3.348, 12.959],
                      "bulk_1T_TiS2": [3.417, 6.419],
                      "2D_2H_TiSe2": [3.400], "2D_1T_TiSe2": [3.544], "bulk_2H_TiSe2": [3.450, 13.000],
                      "bulk_1T_TiSe2": [3.544, 6.694],
                      "2D_2H_TiTe2": [3.560], "2D_1T_TiTe2": [3.779], "bulk_2H_TiTe2": [3.560, 13.200],
                      "bulk_1T_TiTe2": [3.779, 6.817],
                      "2D_2H_VS2": [3.194], "2D_1T_VS2": [3.194], "bulk_2H_VS2": [3.194, 13.000],
                      "bulk_1T_VS2": [3.194, 6.538],
                      "2D_2H_VSe2": [3.354], "2D_1T_VSe2": [3.354], "bulk_2H_VSe2": [3.354, 13.140],
                      "bulk_1T_VSe2": [3.354, 7.000],
                      "2D_2H_VTe2": [3.656], "2D_1T_VTe2": [3.656], "bulk_2H_VTe2": [3.656, 13.400],
                      "bulk_1T_VTe2": [3.656, 6.954],
                      "2D_2H_NbS2": [3.363], "2D_1T_NbS2": [3.377], "bulk_2H_NbS2": [3.363, 13.249],
                      "bulk_1T_NbS2": [3.377, 6.311],
                      "2D_2H_NbSe2": [3.489], "2D_1T_NbSe2": [3.488], "bulk_2H_NbSe2": [3.489, 13.757],
                      "bulk_1T_NbSe2": [3.488, 6.739],
                      "2D_2H_NbTe2": [3.693], "2D_1T_NbTe2": [3.693], "bulk_2H_NbTe2": [3.693, 13.970],
                      "bulk_1T_NbTe2": [3.693, 7.193], "2D_2H_VS2": [3.31, 12.07], "2D_1T_VS2": [3.417, 6.419],
                      "bulk_2H_VS2": [3.31, 12.07], "bulk_1T_VS2": [3.417, 6.419],
                      "2D_2H_TaS2": [3.343], "2D_1T_TaS2": [3.378], "bulk_2H_TaS2": [3.342, 13.760],
                      "bulk_1T_TaS2": [3.378, 6.953],
                      "2D_2H_TaSe2": [3.476], "2D_1T_TaSe2": [3.499], "bulk_2H_TaSe2": [3.476, 14.116],
                      "bulk_1T_TaSe2": [3.499, 6.822],
                      "2D_2H_TaTe2": [3.600], "2D_1T_TaTe2": [3.600], "bulk_2H_TaTe2": [3.600, 14.600],
                      "bulk_1T_TaTe2": [3.600, 6.900],
                      }

reference_cell_TMD = {""}
reference_positions_TMD = {"2H_bulk": [[0.0000000000, 0.0000000000, 0.7500000000],
                                       [0.0000000000, 0.0000000000, 0.2500000000],
                                       [0.6666666667, 0.3333333333, 0.6275000000],
                                       [0.6666666667, 0.3333333333, 0.8725000000],
                                       [0.3333333333, 0.6666666667, 0.1275000000],
                                       [0.3333333333, 0.6666666667, 0.3725000000]
                                       ],
                           "1T_bulk": [[0.0000000000, 0.0000000000, 0.0000000000],
                                       [0.6666666667, 0.3333333333, 0.2200000000],
                                       [0.3333333333, 0.6666666667, 0.7800000000]
                                       ],
                           }

reference_alat_kagome = {"KV3Sb5" :[5.48, 9.37], "RbV3Sb5":[5.50, 9.38],  "CsV3Sb5":[5.49, 9.89],
                         "KNb3Bi5":[5.89, 9.59], "RbNb3Bi5":[5.89, 9.71], "CsNb3Bi5":[5.91, 9.87]}

reference_positions_kagome = [[0.0000000000000000,  0.0000000000000000,  0.0000000000000000],
                              [0.5000000000000000,  0.0000000000000000,  0.5000000000000000],
                              [0.0000000000000000,  0.5000000000000000,  0.5000000000000000],
                              [0.5000000000000000,  0.5000000000000000,  0.5000000000000000],
                              [0.0000000000000000,  0.0000000000000000,  0.5000000000000000],
                              [0.3333333429999996,  0.6666666870000029,  0.7500000000000000],
                              [0.6666666269999979,  0.3333333129999971,  0.2500000000000000],
                              [0.6666666269999979,  0.3333333129999971,  0.7500000000000000],
                              [0.3333333429999996,  0.6666666870000029,  0.2500000000000000]
                             ]


def vc_relax(to_read_orig_poscar):
    ifile = 0

    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:

            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("02_vc_relax"):
                print("VC relaxation file exists: %s" % ("02_vc_relax"))
            else:
                print("Creating post VC file for bulk 2H %s%s2 \n" % (cation, anion))
                os.mkdir("02_vc_relax")

            os.chdir("./02_vc_relax")
            os.system("cp ../02_volume/POTCAR ./")
            os.system("cp ../../%s/POSCAR_%s%s2 ./POSCAR" % (to_read_orig_poscar, cation, anion) )
            os.system("cp ../../../inputs/INCAR_vc_relax ./INCAR")
            os.system("cp ../../../inputs/KPOINTS_relax ./KPOINTS")
            os.system("cp ../../../inputs/vc_relax.sh ./")

            os.system("sbatch --gpus=1 vc_relax.sh")
            print(" %s%s2 done !\n" % (cation, anion))

            ifile += 1

def relax(to_read_orig_poscar):
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:


            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("03_relax"):
                print("Post relaxation file exists: %s" % ("03_relax"))
            else:
                print("Creating post relaxation file for bulk 2H %s%s2 \n" % (cation, anion))
                os.system("mkdir 03_relax")

            os.chdir("./03_relax")
            os.system("cp ../02_vc_relax/POTCAR ./")
            os.system("cp ../../%s/POSCAR_%s%s2 ./POSCAR" % (to_read_orig_poscar, cation, anion) )
            os.system("cp ../../../inputs/INCAR_relax ./INCAR")
            os.system("cp ../../../inputs/KPOINTS_relax ./KPOINTS")
            os.system("cp ../../../inputs/relax.sh ./")
            os.system("sbatch --gpus=1 relax.sh")

            ifile += 1

def check_relaxation(outdir, outrelax):
    ifile = 0
    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion) + "%s/" %outdir)
            os.system("cp ./CONTCAR ../../relaxed/POSCAR_%s%s2" %(cation, anion))
            if os.path.isfile(outrelax):
                print("Open file: %s \n" %(outrelax))
                with open(outrelax) as f:
                    lines = f.readlines()
                    for line in lines:
                        if "reached required accuracy" in line:
                            print("%02d_%s%s2 reached accuracy !\n" % (ifile, cation, anion))

                with open("OUTCAR") as f:
                    lines = f.readlines()
                    stress = []
                    for iline in range(len(lines)):
                        if "FORCE on cell" in lines[iline]:
                            stress.append(lines[iline + 13])
                    print(stress[-1])
            else:
                print("%s NOT FOUND ! \n" %( outrelax ))
            ifile += 1


def scf():
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("04_scf"):
                print("SCF file exists: %s" % ("04_scf"))
            else:
                os.system("mkdir 04_scf")

            os.chdir("./04_scf")
            os.system("cp ../03_relax/POTCAR ./")
            os.system("cp ../../relaxed/POSCAR_%s%s2 ./POSCAR" %(cation, anion))
            os.system("cp ../../../inputs/INCAR_scf ./INCAR")
            os.system("cp ../../../inputs/KPOINTS_scf ./KPOINTS")
            os.system("cp ../../../inputs/scf.sh ./")

            if os.path.isfile("POSCAR") and os.path.isfile("POTCAR") and os.path.isfile("INCAR") and os.path.isfile("KPOINTS"):
                print("Inputs file ready ! \n")
            else:
                print("Inputs file NOT complete ! \n")
            os.system("sbatch --gpus=1 scf.sh")

            os.chdir(root)
            ifile += 1

def check_scf(outscf):
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))

            if os.path.isdir(root + "/%02d_%s%s2/" % (ifile, cation, anion) + "04_scf"):
                print("SCF file exists: %s \n" % ("04_scf"))
            else:
                print("SCF NOT FOUND for bulk 2H %s%s2 \n" % (cation, anion))
                break
            if os.path.isfile(root + "/%02d_%s%s2/" % (ifile, cation, anion) + "04_scf/" + outscf):
                print("Open file: %s \n" % (outscf))
                with open(root + "/%02d_%s%s2/" % (ifile, cation, anion) + "04_scf/" + outscf) as f:
                    lines = f.readlines()
                    for line in lines:
                        if "F=" in line:
                            print("%02d_%s%s2 SCF done !\n" % (ifile, cation, anion))

            os.chdir(root)
            ifile += 1


def dos():
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:

            with open(root + "/%02d_%s%s2/04_scf/DOSCAR" % (ifile, cation, anion), "r") as f:
                lines = f.readlines()
                FERMI_ENERGY = float(lines[5].strip().split()[3])

            os.chdir(root + "/%02d_%s%s2" % (ifile, cation, anion))
            if os.path.isdir("08_dos"):
                print("SCF file exists: %s" % ("08_dos"))
            else:
                os.system("mkdir 08_dos")

            os.chdir("./08_dos")
            os.system("cp ../04_scf/POSCAR ./")
            os.system("cp ../04_scf/POTCAR ./")
            os.system("cp ../04_scf/CHG ./")
            os.system("cp ../04_scf/CHGCAR ./")

            os.system("cp ../../../inputs/INCAR_scf ./INCAR")
            os.system("cp ../../../inputs/KPOINTS_scf ./KPOINTS")
            os.system("cp ../../../inputs/scf.sh ./")

            with open("INCAR", "a") as f:
                f.write("\n")
                f.write("NEDOS   =  1000 \n")
                f.write("EMIN    =  %.2f \n" %(FERMI_ENERGY - 2.0))
                f.write("EMAX    =  %.2f \n" %(FERMI_ENERGY + 2.0))


            if os.path.isfile("POSCAR") and os.path.isfile("POTCAR") and os.path.isfile("INCAR") and os.path.isfile("KPOINTS"):
                print("Inputs file ready ! \n")
            else:
                print("Inputs file NOT complete ! \n")
            os.system("sbatch --gpus=1 scf.sh")

            ifile += 1

def collect_dos():
    root = os.getcwd()
    ifile = 0
    if os.path.isdir(root + "/DOS"):
        print("FOUND DoS file: %s\n" % (root + "/DOS"))
    else:
        print("DoS file NOT FOUND \n Creating new file: %s \n" %(root + "/DOS"))
        os.mkdir(root + "/DOS")

    for cation in cations:
        for anion in anions:
            print("Collecting DOSCAR for %s%s2\n" % (cation, anion))
            if os.path.isdir(root + "/DOS" + "/%02d_%s%s2" % (ifile, cation, anion)):
                print("FOUND DoS file for /%02d_%s%s2\n" % (ifile, cation, anion))
            else:
                os.mkdir(root + "/DOS" + "/%02d_%s%s2" % (ifile, cation, anion))

            os.chdir(root + "/%02d_%s%s2/08_dos" % (ifile, cation, anion))
            os.system("cp DOSCAR %s/DOS/%02d_%s%s2" %(root, ifile, cation, anion))
            print("%s%s2 Done!\n" % (cation, anion))

            ifile += 1


def band():
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("05_band"):
                print("Band structure file exists: %s" % ("05_band"))
            else:
                print("Creating band structure file for: %02d_%s%s2" % (ifile, cation, anion))
                os.system("mkdir 05_band")

            os.chdir("./05_band")
            os.system("cp ../04_scf/POSCAR ./")
            os.system("cp ../04_scf/POTCAR ./")
            os.system("cp ../04_scf/CHG ./")
            os.system("cp ../04_scf/CHGCAR ./")
            os.system("cp ../../../inputs/INCAR_band ./INCAR")
            os.system("cp ../../../inputs/KPOINTS_band ./KPOINTS")
            os.system("cp ../../../inputs/band.sh ./")
            os.system("sbatch --gpus=1 band.sh")

            print("%02d_%s%s2 band done !\n" % (ifile, cation, anion))
            os.chdir(root)
            ifile += 1


def collect_band_for_vasp():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/band_for_vaspkit"):
        print("Band data file exists: band_for_vaspkit \n")
    else:
        print("Creating band data file \n")
        os.mkdir(root + "/band_for_vaspkit")

    for cation in cations:
        for anion in anions:

            if os.path.isdir(root + "/band_for_vaspkit/%02d_%s%s2/" % (ifile, cation, anion)):
                print("Found band file for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating band file for: %02d_%s%s2 \n" % (ifile, cation, anion))
                os.mkdir(root + "/band_for_vaspkit/%02d_%s%s2" % (ifile, cation, anion))
            os.chdir(root + "/%02d_%s%s2/04_scf/" % (ifile, cation, anion))
            os.system("cp DOSCAR %s" % (
                        root + "/band_for_vaspkit/%02d_%s%s2/" % (ifile, cation, anion)))
            os.chdir(root + "/%02d_%s%s2/05_band/" % (ifile, cation, anion))
            os.system("cp -r INCAR POSCAR KPOINTS POTCAR EIGENVAL PROCAR %s" %(root + "/band_for_vaspkit/%02d_%s%s2/"  % (ifile, cation, anion)))

            ifile += 1

def collect_band_dat():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/band_dat"):
        print("Band data file exists: band_dat \n")
    else:
        print("Creating band data file \n")
        os.mkdir(root + "/band_dat")

    for cation in cations:
        for anion in anions:
            if os.path.isdir(root + "/band_dat/%02d_%s%s2" % (ifile, cation, anion)):
                print("Band data file exists for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating band data file for: %02d_%s%s2  \n" % (ifile, cation, anion))
                os.mkdir(root + "/band_dat/%02d_%s%s2" % (ifile, cation, anion))


            os.chdir(root + "/band_for_vaspkit/%02d_%s%s2" % (ifile, cation, anion))

            os.system("cp -r FERMI_ENERGY KLABELS KLINES.dat PBAND_%s.dat PBAND_%s.dat ../../band_dat/%02d_%s%s2" % (cation, anion, ifile, cation, anion))

            ifile += 1


def fs_ibz(kptfile):
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("06_fs_ibz"):
                print("Fermi surface for irreducible Brillouin zone file exists: %s" % ("06_fs_ibz"))
            else:
                os.system("mkdir 06_fs_ibz")

            os.chdir("./06_fs_ibz")
            os.system("cp ../04_scf/POSCAR ./")
            os.system("cp ../04_scf/POTCAR ./")
            os.system("cp ../04_scf/CHG ./")
            os.system("cp ../04_scf/CHGCAR ./")
            os.system("cp ../../../inputs/INCAR_fs ./INCAR")
            os.system("cp ../../../inputs/%s ./KPOINTS" %kptfile)
            os.system("cp ../../../inputs/fs.sh ./")
            os.system("sbatch --gpus=1 fs.sh")
            os.chdir(root)
            ifile += 1



def collect_fs_ibz():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/fs_ibz"):
        print("Fermi surface data file exists: fs_ibz \n")
    else:
        print("Creating Fermi surface data file \n")
        os.mkdir(root + "/fs_ibz")


    for cation in cations:
        for anion in anions:
            if os.path.isdir(root + "/fs_ibz/%02d_%s%s2" % (ifile, cation, anion)):
                print("Fermi surface data file exists for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating Fermi surface data file for: %02d_%s%s2 \n" % (ifile, cation, anion))
                os.mkdir(root + "/fs_ibz/%02d_%s%s2" % (ifile, cation, anion))
            os.chdir(root + "/%02d_%s%s2/04_scf" % (ifile, cation, anion))
            os.system("cp  DOSCAR %s/fs_ibz/%02d_%s%s2/" % (root, ifile, cation, anion))

            os.chdir(root + "/%02d_%s%s2/06_fs_ibz" % (ifile, cation, anion))

            os.system("cp -r EIGENVAL KPOINTS %s/fs_ibz/%02d_%s%s2/" %(root, ifile, cation, anion))

            print("%02d_%s%s2/ Done ! \n" % (ifile, cation, anion))

            ifile += 1


def fs_fbz(kptfile):
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("07_fs_fbz"):
                print("Fermi surface for full Brillouin zone file exists: %s" % ("07_fs_fbz"))
            else:
                print("Creating file for FBZ Fermi surface: %s" % ("07_fs_fbz"))
                os.system("mkdir 07_fs_fbz")

            os.chdir("./07_fs_fbz")
            os.system("cp ../04_scf/POSCAR ./")
            os.system("cp ../04_scf/POTCAR ./")
            os.system("cp ../04_scf/CHG ./")
            os.system("cp ../04_scf/CHGCAR ./")
            os.system("cp ../../../inputs/INCAR_fs ./INCAR")
            os.system("cp ../../../inputs/%s ./KPOINTS" %(kptfile))
            os.system("cp ../../../inputs/fs_fbz.sh ./")
            # os.system("sbatch --gpus=2 fs_fbz.sh")

            print("%02d_%s%s2/ Done ! \n" % (ifile, cation, anion))

            os.chdir(root)
            ifile += 1

def collect_fs_fbz():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/fs_fbz"):
        print("Fermi surface data file exists: fs_fbz \n")
    else:
        print("Creating Fermi surface data file \n")
        os.mkdir(root + "/fs_fbz")


    for cation in cations:
        for anion in anions:
            if os.path.isdir(root + "/fs_fbz/%02d_%s%s2" % (ifile, cation, anion)):
                print("Fermi surface data file exists for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating Fermi surface data file for: %02d_%s%s2 \n" % (ifile, cation, anion))
                os.mkdir(root + "/fs_fbz/%02d_%s%s2" % (ifile, cation, anion))
            os.chdir(root + "/%02d_%s%s2/04_scf" % (ifile, cation, anion))
            os.system("cp  DOSCAR %s/fs_fbz/%02d_%s%s2/" % (root, ifile, cation, anion))

            os.chdir(root + "/%02d_%s%s2/07_fs_fbz" % (ifile, cation, anion))

            os.system("cp -r EIGENVAL DOSCAR KPOINTS POSCAR %s/fs_fbz/%02d_%s%s2/" %(root, ifile, cation, anion))

            print("%02d_%s%s2/ Done ! \n" % (ifile, cation, anion))

            ifile += 1


def collect_poscar_for_phonopy():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/poscar_for_phonopy"):
        print("POSCAR for phonopy file exists: poscar_for_phonopy \n")
    else:
        print("Creating POSCAR for phonopy file \n")
        os.mkdir(root + "/poscar_for_phonopy")

    for cation in cations:
        for anion in anions:
            if os.path.isdir(root + "/poscar_for_phonopy/%02d_%s%s2/" % (ifile, cation, anion)):
                print("Found band file for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating POSCAR for phonopy file for: %02d_%s%s2 \n" % (ifile, cation, anion))
                os.mkdir(root + "/poscar_for_phonopy/%02d_%s%s2" % (ifile, cation, anion))

            os.chdir(root + "/%02d_%s%s2/04_scf/" % (ifile, cation, anion))
            os.system("cp POSCAR  %s/poscar_for_phonopy/%02d_%s%s2/" %(root, ifile, cation, anion))

            ifile += 1


def phonon():
    ifile = 0
    print("Presen file: %s" % root)

    for cation in cations:
        for anion in anions:
            os.chdir(root + "/%02d_%s%s2/" % (ifile, cation, anion))
            if os.path.isdir("09_phonon"):
                print("phonon file exists: 09_phonon")
            else:
                os.mkdir("./09_phonon")


            os.system("cp ./04_scf/POTCAR ./09_phonon")
            os.system("cp ../../inputs/INCAR_phonon ./09_phonon/INCAR")
            os.system("cp ../../inputs/KPOINTS_phonon ./09_phonon/KPOINTS")
            os.system("cp ../../inputs/phonon.sh ./09_phonon/")
            if os.path.isfile("../poscar_for_phonopy/%02d_%s%s2/SPOSCAR" % (ifile, cation, anion)):
                print("Found super poscar for %s%s2 \n" % (cation, anion))
                os.system("cp ../poscar_for_phonopy/%02d_%s%s2/SPOSCAR ./09_phonon/POSCAR" % (ifile, cation, anion))
            else:
                print("No super POSCAR !\n")
                ifile += 1
                continue
            os.chdir("./09_phonon")
            if os.path.isfile("POSCAR") and os.path.isfile("POTCAR") and os.path.isfile("INCAR") and os.path.isfile("KPOINTS"):
                print("Inputs file ready ! \n")
            else:
                print("Inputs file NOT complete ! \n")
            os.system("sbatch --gpus=1 phonon.sh")

            ifile += 1


def collect_phonon_xml():
    ifile = 0
    print("Presen file: %s" % root)

    if os.path.isdir(root + "/xml_for_phonon"):
        print("xml file exists: xml_for_phonon \n")
    else:
        print("Creating xml_for_phonon file \n")
        os.mkdir(root + "/xml_for_phonon")

    for cation in cations:
        for anion in anions:
            if os.path.isdir(root + "/xml_for_phonon/%02d_%s%s2/" % (ifile, cation, anion)):
                print("Found band file for: %02d_%s%s2 \n" % (ifile, cation, anion))
            else:
                print("Creating POSCAR for phonopy file for: %02d_%s%s2 \n" % (ifile, cation, anion))
                os.mkdir(root + "/xml_for_phonon/%02d_%s%s2" % (ifile, cation, anion))

            os.chdir(root + "/%02d_%s%s2/09_phonon/" % (ifile, cation, anion))
            os.system("cp -r ../04_scf/POSCAR vasprun.xml %s/xml_for_phonon/%02d_%s%s2/" % (root, ifile, cation, anion))

            ifile += 1
# cleanall()
# vc_relax("04_scf")
# scf()
# check_scf("outscf")
# band()
# check_relaxation("03_relax", "outrelax")
# collect_band_for_vasp()
# collect_fs_ibz()
# fs_fbz("KPOINTS_45_45_1")
# collect_band_dat()
# collect_fs_fbz()

# collect_dos()
# dos()
# fs_ibz("KPOINTS_ibz_zhalf")
# collect_poscar_for_phonopy()
collect_phonon_xml()