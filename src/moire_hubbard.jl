using ITensors
using ITensorMPS
using HDF5
using Plots
using LaTeXStrings
using LinearAlgebra
using Random

# Disable BLAS and Strided multithreading, as it conflicts with Julia threads.
BLAS.set_num_threads(1)
ITensors.Strided.disable_threads()

ITensors.disable_debug_checks()

"""
Finite quasi-1D moiré Hubbard model with open boundary conditions.

INPUTS: L is the y-axis length of the ladder; nsweeps is the number of DMRG sweeps;
maxdim is the maximum bond dimension kept after each sweep, ie. generally you
want to increase the number kept as you do more sweeps; cutoff is the sought-after
truncation error in each sweep (may be array for each sweep); t is the 
hopping constant; U_0 is the on-site density-density interaction coupling;
U_1 is the nearest-neighbour density-density interaction coupling; mu is the chemical
potential; phi is (kind-of) a magnetic field.

RETURNS: energy eigenvalues of psi0, truncated ground states psi0.
"""
function moire_hubbard(L, nsweeps, maxdim, cutoff, t, U_0, U_1, mu, phi)

    N = 2L
    hop = -t*exp(-im*phi)
    hop_conj = -t*exp(im*phi)

    # Pre-allocate result containers (thread-safe: each thread writes to
    # a distinct index; push! would be a data race on the array length field).
    energies = Vector{Float64}(undef, length(mu))
    groundstates = Vector{MPS}(undef, length(mu))

    # Create N physical indices w/ the electron property tag.
    sites = siteinds("Electron", N)

    print("Building Hamiltonian...")
    # Build the mu-independent part of the Hamiltonian.
    os0 = OpSum()

    # Nearest neighbour terms.
    for j in 1:(N-1)
        os0 += hop, "Cdagup", j, "Cup", j+1
        os0 += hop_conj, "Cdagup", j+1, "Cup", j
        os0 += hop_conj, "Cdagdn", j, "Cdn", j+1
        os0 += hop, "Cdagdn", j+1, "Cdn", j
    end

    # Second nearest neighbour terms.
    for j in 1:(N-2)
        os0 += hop_conj, "Cdagup", j, "Cup", j+2
        os0 += hop, "Cdagup", j+2, "Cup", j
        os0 += hop, "Cdagdn", j, "Cdn", j+2
        os0 += hop_conj, "Cdagdn", j+2, "Cdn", j
    end

    # On-site density-density interactions
    for j in 1:N
        os0 += U_0, "Nup", j, "Ndn", j
    end

    # Nearest-neighbour density-density interactions.
    for j in 1:(N-1)
        os0 += U_1, "Ntot", j, "Ntot", j+1
    end
    for j in 1:(N-2)
        os0 += U_1, "Ntot", j, "Ntot", j+2
    end

    H0 = MPO(os0, sites)

    # Chemical potential term.
    os_mu = OpSum()
    for n in 1:N
        os_mu += -1.0, "Ntot", n
    end
    H_mu = MPO(os_mu, sites)
    
    print("\nBeginning DMRG...")
    # Threads.@threads over the independent mu values. Each task gets its own
    # fresh randomMPS initial guess (no shared mutable state) and writes to
    # distinct pre-allocated slots (no data race). Run with julia --threads=N.
    Threads.@threads for i in eachindex(mu)
        m = mu[i]

        print("\n------------------ mu = $m ------------------\n")

        # Seeded per-index RNG so results are reproducible run-to-run.
        Random.seed!(i)

        # Create a random initial guess for the ground state MPS
        # with sites as physical indices and with bond dimension 10. 
        psi0 = randomMPS(sites, 10) 

        # Assemble the full Hamiltonian.
        H = H0 + m*H_mu

        energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
        energies[i] = energy
        groundstates[i] = psi

    end

    return energies, groundstates

end