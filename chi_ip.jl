using LinearAlgebra
using Printf
using DataFrames
using CSV
using Plots

current_file = pwd()

# cations = ["Nb", "Ta"]
# anions = ["S", "Se", "Te"]

cations = ["Ti", "V", "Nb", "Ta"]
anions = ["S", "Se", "Te"]

Hartree2eV = 27.211386245988
Boltzmann_eV = 8.617333262145 * 1e-5 # eV
Bohr_radius = 0.52917721067 # ang

Nkx = Nky = 63
Nqx = Nqy = 63
Nqz = 1
Nkz = 1
RECIPROCAL_LATTICE = [sqrt(3)/2    1/2          0.000000000;
                      0.000000000  1.000000000  0.000000000;
                      0.000000000  0.000000000  0.028097780]

SYMMETRY_ops = cat([ 1  0  0;  0  1  0; 0  0  1],
                   [ 0  1  0;  1  0  0; 0  0  1],
                   [ 1  1  0;  0 -1  0; 0  0  1],
                   [ 0 -1  0;  1  1  0; 0  0  1],
                   [ 1  1  0; -1  0  0; 0  0  1],
                   [ 1  0  0; -1 -1  0; 0  0  1],
                   [-1  0  0;  0 -1  0; 0  0  1],
                   [ 0 -1  0; -1  0  0; 0  0  1],
                   [-1 -1  0;  0  1  0; 0  0  1],
                   [ 0  1  0; -1 -1  0; 0  0  1],
                   [-1 -1  0;  1  0  0; 0  0  1],
                   [-1  0  0;  1  1  0; 0  0  1], dims=3)




function Fermi_Dirac(E, T=300)
    return 1 / (1 + exp(E/ (Boltzmann_eV * T)))
end

function delta1(x, ϵ)
    return ϵ / (pi * (x ^ 2 + ϵ ^ 2))
end

function delta1(x, ϵ)
    return exp(- x^2 / ϵ) / (2 * sqrt(pi * ϵ))
end

function delta3(x, ϵ)
    return ϵ * (abs(x) ^ (ϵ - 1))
end


function load_kpoints(filname)
    KPOINTS_frac = []
    NKPOINTS = 0
    @printf "Loading kpoints from file: %s \n" filname
    open(filname, "r") do kptio
        lines = readlines(kptio)
        NKPOINTS = parse(Int32, strip(lines[2]))
        for line in lines[4:length(lines)]
            line = [parse(Float64, ss) for ss in split(line)]
            append!(KPOINTS_frac, line)
        end
    end

    KPOINTS_frac = transpose(reshape(KPOINTS_frac, (4,NKPOINTS)))
    return NKPOINTS, KPOINTS_frac
end

function load_fermi_energy(filname)
    FERMI_ENERGY = 0.0
    @printf "Loading Fermi energies from file: %s \n" string(filname, "\\DOSCAR")
    open(string(filname), "r") do fermiio
        line = split(readlines(fermiio)[6])
        FERMI_ENERGY = parse(Float64, line[4])
    end 

    @printf "Fermi energy: %.4f ! \n" FERMI_ENERGY
    return FERMI_ENERGY
end

function load_energy(filname, FERMI_ENERGY)
    NELEC, NKPOINTS, NBANDS = 0, 0, 0

    ENERGIES = []
    KPOINTS_frac = []

    @printf "Loading electron energies from file: %s \n" string(filname, "\\EIGENVAL")
    open(string(filname, "\\EIGENVAL"), "r") do energyio
        lines = readlines(energyio)
        NELEC, NKPOINTS, NBANDS = [parse(Int32, ss) for ss in split(lines[6])]
        for line in lines[7:length(lines)]
            line = [parse(Float64, ss) for ss in split(line)]
            if length(line) == 4
                append!(KPOINTS_frac, line)
            end
            if length(line) == 3
                append!(ENERGIES, line)
            end
        end
    end 

    ENERGIES = transpose(reshape(ENERGIES, (3, NBANDS, NKPOINTS))[2, :, :]) .- FERMI_ENERGY
    KPOINTS_frac = transpose(reshape(KPOINTS_frac, (4, NKPOINTS)))

    @printf "Energy done ! \n"
    return NKPOINTS, NBANDS, ENERGIES, KPOINTS_frac

end

function find_fermi_bands(ENERGIES, NBANDS)
    FERMI_BAND = []
    for iband in 1:NBANDS
        lower_bound, upper_bound = minimum(ENERGIES[:, iband]), maximum(ENERGIES[:, iband])

        if (lower_bound < 0) && (upper_bound > 0)
            @printf "Metal band %d ranging from %.4f to %.4f \n" iband lower_bound upper_bound
            append!(FERMI_BAND, [iband, lower_bound, upper_bound])
        end
    end
    FERMI_BAND = reshape(FERMI_BAND, (3, length(FERMI_BAND) ÷ 3))
    return FERMI_BAND
end


function calculate_chi_ip(NQPOINTS, QPOINTS, NKPOINTS, KPOINTS, ENERGIES, FERMI_BAND, T=10, ϵ=1e-3)
    chi_ip = zeros(Float64, 2, size(FERMI_BAND)[2], NQPOINTS)
    for iq in 1:NQPOINTS
        begin_time = time()
        for iband in FERMI_BAND[1,:]
            for ik in 1:NKPOINTS
                for jband in FERMI_BAND[1,:]
                    kq = (QPOINTS[iq, 1:3] + KPOINTS[ik, 1:3] + [1, 1, 1]) .% 1
                    kqid = 1 + Int(round(kq[1] * Nkx)) % Nkx + Nkx * (Int(round(kq[2] * Nky)) % Nky)
                    Ek = ENERGIES[ik, Int(iband)]
                    Ekq = ENERGIES[kqid, Int(jband)]

                    if abs(Ekq - Ek) < 1e-8
                        Ekq += 1e-8
                    end
                    if abs(Ekq) < 1e-8
                        Ekq += 1e-8
                    end
                    if abs(Ek) < 1e-8
                        Ek += 1e-8
                    end
                    fk  = Fermi_Dirac(Ek, T)
                    fkq = Fermi_Dirac(Ekq, T)
                    chi_ip[1, Int(iband - FERMI_BAND[1,1] + 1), iq] += (fkq - fk) / ((Ekq - Ek) * Nkx * Nky)
                    chi_ip[2, Int(iband - FERMI_BAND[1,1] + 1), iq] +=  delta1(Ek, ϵ) * delta1(Ekq, ϵ)
                end
            end
        end
        consume = time() - begin_time
        remaining = (NQPOINTS - iq) * consume
        @printf "%d / %d finished in %.4f s! Approximately %.4f s remaining \n" iq NQPOINTS consume remaining
    end
    chi_ip = chi_ip ./ (Nkx * Nky * Nkz)
    return chi_ip
end


if isdir("chi_ip")
    @printf "Found χ directory: %s" string(current_file, "\\chi_ip")
else
    mkdir(string(current_file, "\\chi_ip"))
end



global imat = 0
for cation in cations
    for anion in anions
        filfermi = @sprintf "%s\\DOS\\%02d_%s%s2\\DOSCAR" current_file imat cation anion
        FERMI_ENERGY = load_fermi_energy(filfermi)
        filkpts_ibz = @sprintf "%s\\chi_ip\\KPOINTS_63_63_1" current_file
        NQPOINTS, QPOINTS_frac = load_kpoints(filkpts_ibz)
        filfs_fbz = @sprintf "%s\\fs_fbz\\%02d_%s%s2" current_file imat cation anion
        NKPOINTS, NBANDS, ENERGIES, KPOINTS_frac = load_energy(filfs_fbz, FERMI_ENERGY)

        FERMI_BAND = find_fermi_bands(ENERGIES, NBANDS)

        chi_ip = calculate_chi_ip(NQPOINTS, QPOINTS_frac, NKPOINTS, KPOINTS_frac, ENERGIES, FERMI_BAND)

        
        chi_ip_file = @sprintf "%s\\chi_ip\\%02d_%s%s2" current_file imat cation anion 


        open(chi_ip_file, "w") do chireio
            for iband in FERMI_BAND[1,:]
                @printf(chireio, "Reχ(%d)       Imχ(%d)       ", iband, iband)
            end

            @printf(chireio, "\n")
            for iq in 1:NQPOINTS
                for iband in FERMI_BAND[1,:]
                    @printf(chireio, "%.8f    %.8f    ", chi_ip[1, Int(iband - FERMI_BAND[1,1] + 1), iq], chi_ip[2, Int(iband - FERMI_BAND[1,1] + 1), iq])
                end
                @printf(chireio, "\n")
            end
        end

        global imat += 1
    end
end
