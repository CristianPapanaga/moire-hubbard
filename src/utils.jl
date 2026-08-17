using ITensors
using ITensorMPS
using HDF5
using Plots; pgfplotsx()
default(legend_font_halign=:left) # By default legend entries are centered; left-align them.
using Plots.Measures
using LaTeXStrings
using NaturalSort
using Base.Iterators
using FFTW

"""
Given lattice vectors a, b, generates a N_x * N_y lattice of unit cells, where
the sites are at coordinates site_unit_cell_coords within the unit cell, given in
terms of lattice vectors a, b.
"""
function lattice_coords_2d(a, b, site_unit_cell_coords, N_x, N_y)

    lattice_x_coords = []
    lattice_y_coords = []

    for site in site_unit_cell_coords
        # Define initial x and y unit cell coordinates.
        x_initial = site[1]*a[1] + site[2]*b[1]
        y_initial = site[2]*a[2] + site[2]*b[2]

        # Generate unit cell coordinates for each pair of points indexing the site,
        # as given by the Cartesian product.
        x_coords = [(x_initial + i*a[1] + j*b[1]) for (i, j) in Base.product(1:N_x, 1:N_y)]
        y_coords = [(y_initial + i*a[2] + j*b[2]) for (i, j) in Base.product(1:N_x, 1:N_y)]

        push!(lattice_x_coords, x_coords)
        push!(lattice_y_coords, y_coords)
    end

    return (collect(flatten(lattice_x_coords')), collect(flatten(lattice_y_coords')))
    
end

"""
Load the simulation parameters and resulting calculated ground states
from file_name.

RETURNS: ground_states, L, t, U_0, U_1, mu, phi.
"""
function load_data(file_name)

    # Access h5 file containing MPS data.
    file = h5open(file_name, "r")

    # Retrieve names of datasets corresponding to MPS
    # and store list of ground states.
    mps_names = sort(keys(file), lt = natural)
    ground_states = []
    for mps_name in mps_names

        # I open and close the file every single
        # iteration because there seems to be some bug
        # in the HDF5 C library where each consecutive
        # read operation takes longer and longer; closing
        # the file resets it, I guess? Idk it works.

        file = h5open(file_name, "r")

        psi0 = read(file, "$mps_name", MPS)
        push!(ground_states, psi0)

        close(file)

    end

    # Retrieve system parameters from file name & MPS dataset name.
    regex = r"(?<=\[)(.*?)(?=\])" # All values between [ and ], excluding the delimiters.
    L, t, U_0, U_1, phi = [parse(Float64, i.match) for i in eachmatch(regex, file_name)]
    mu = []
    for mps_name in mps_names
        m = [parse(Float64, i.match) for i in eachmatch(regex, mps_name)]
        push!(mu, m[1])
    end  

    return ground_states, L, t, U_0, U_1, mu, phi

end

"""
Plot density as a function of chemical potential for a single dataset.
"""
function plot_density_vs_mu(ground_states, L, t, U_0, U_1, mu, phi)

    densities = []
    for psi0 in ground_states

        Ntot = expect(psi0, "Ntot")
        number = sum(Ntot)
        density = number/2L

        push!(densities, density)

    end

    plot(mu, densities, lc = :red, minorgrid = true,
     label = latexstring("U = $U_0"), legend = :topleft)
    scatter!(mu, densities, mc = :black, markersize = 1, label = "DMRG")
    xlabel!(L"\mu")
    ylabel!(L"\rho")

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")
    savefig("$folder/Density L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi].pdf")

end

"""
Plot density as a function of chemical potential for all values of U
in the given folder.
"""
function plot_density_vs_mu_all_U(folder, L, t, phi)

    # Initialize plot.
    plot(legend = :topleft, minorgrid = true, color_palette = :tab10) 
    xlabel!(L"\mu")
    ylabel!(L"\rho")

    # Get all file names in folder
    # and plot the density vs mu.
    file_names = readdir(folder, sort = false)
    file_names = sort(file_names, lt = natural)

    for file_name in file_names
    
        ground_states, L, t, U_0, U_1, mu, phi = load_data("$folder/$file_name")

        densities = []
        for psi0 in ground_states

            Ntot = expect(psi0, "Ntot")
            number = sum(Ntot)
            density = number/2L

            push!(densities, density)

        end

        plot!(mu, densities, label = latexstring("U = $U_0"))

    end

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")
    savefig("$folder/Density for different U L = [$L]  t = [$t]  phi = [$phi].pdf")

end

"""
Plot the real space lattice occupation for a given mu index.
"""
function plot_lattice_occupation(ground_states, L, t, U_0, U_1, mu, phi, mu_index, x, y)

    current_mu = mu[mu_index]

    # Colour lattice sites based on occupation expectation.
    Ntots = []
    for psi0 in ground_states

        Ntot = expect(psi0, "Ntot")
        push!(Ntots, Ntot)

    end

    scatter(x, y, marker_z = Ntots[mu_index], color = reverse(cgrad(:devon)),
            legend = false, colorbar = :top, markerstrokewidth = 0,
            clims = (0, 2),
            colorbar_title = L"Total occupation number $\left<n_i^\dagger n_i\right>$",
            aspect_ratio = 3, grid = false, xticks = false, yticks = false, axis = false,
            size = (600, 140))

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")  
    savefig("$folder/Real space occupation L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1] mu = [$current_mu]  phi = [$phi].pdf")

end

"""
Plot the S_z spin correlations for a given mu index.
"""
function plot_spin_correlations(ground_states, L, t, U_0, U_1, mu, phi, mu_index, x, y)

    current_mu = mu[mu_index]

    # Colour lattice sites based on spin correlation.
    spin_corr_matrices = []
    for psi0 in ground_states

        # We take the real part since it has imaginary part
        # which I assume is an error.
        spin_corr_matrix = real.(correlation_matrix(psi0, "Sz", "Sz"))
        push!(spin_corr_matrices, spin_corr_matrix)

    end

    # Pick a central site (site L) as the spin correlation reference
    # and extract the correlation values of site L and every other site
    # in the system.
    spin_corrs = spin_corr_matrices[mu_index][Int(L), 1:end .!= L]

    # When the spin_corrs are too small, plots writes tick labels in scientific
    # notation, which is very ugly; pgfplotsx doesn't currently work with plots' 
    # formatter, so to prevent this scale the spin_corrs and plot those, stating the 
    # scale factor in the title.
    scale_exponent = 2 
    spin_corrs_scaled = spin_corrs .* 10^(scale_exponent)
    scatter(x[1:end .!= L], y[1:end .!= L], marker_z = spin_corrs_scaled,
            color = cgrad([:blue, :white]),
            legend = false, colorbar = :top, markerstrokewidth = 0,
            colorbar_title = latexstring("Spin correlation \$\\left<S^z_i S^z_L\\right>\$ (\$\\times 10^{$scale_exponent}\$)"),
            aspect_ratio = 3, grid = false, xticks = false, yticks = false, axis = false,
            size = (600, 140))

    # Draw cross at reference point.
    scatter!([x[Int(L)]], [y[Int(L)]], markershape = :xcross, markersize = 6, markerstrokewidth = 4)

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")
    savefig("$folder/Spin correlations L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1] mu = [$current_mu]  phi = [$phi].pdf")

end

"""
Plot the charge structure factor (the 1d Fourier transform of Ntot) for a given mu index.
"""
function plot_charge_structure_factor(ground_states, L, t, U_0, U_1, mu, phi, mu_index)

    current_mu = mu[mu_index]

    # Charge structure factor - 1d FT of Ntot.
    Ntots = []
    for psi0 in ground_states

        Ntot = expect(psi0, "Ntot")
        push!(Ntots, Ntot)

    end

    S_c = 1/(2L) .* fftshift(fft(Ntots[mu_index]))
    k = fftshift(fftfreq(Int(2L), Int(2L)))

    plot(k, abs.(S_c), lc = :red, minorgrid = true, legend = false)
    scatter!(k, abs.(S_c), mc = :black, markersize = 2)
    xlabel!(L"k_x")
    ylabel!(L"S_c\left(k_x\right)")

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")
    savefig("$folder/Charge structure factor L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1] mu = [$current_mu]  phi = [$phi].pdf")

end

"""
Plot the magnetic structure factor (the 1d Fourier transform of spin_corr) for a given mu index.
"""
function plot_magnetic_structure_factor(ground_states, L, t, U_0, U_1, mu, phi, mu_index)

    current_mu = mu[mu_index]

    # Magnetic structure factor - 1d FT of spin_corr.
    spin_corr_matrices = []
    for psi0 in ground_states

        # We take the real part since it has imaginary part
        # which I assume is an error.
        spin_corr_matrix = real.(correlation_matrix(psi0, "Sz", "Sz"))
        push!(spin_corr_matrices, spin_corr_matrix)

    end

    spin_corrs = spin_corr_matrices[mu_index][Int(L), 1:end .!= L]

    S_m = 1/(2L) .* fftshift(fft(spin_corrs))
    k = fftshift(fftfreq(Int(2L-1), Int(2L-1)))

    plot(k, abs.(S_m), lc = :red, minorgrid = true, legend = false)
    scatter!(k, abs.(S_m), mc = :black, markersize = 2)
    xlabel!(L"k_x")
    ylabel!(L"S_m\left(k_x\right)")

    folder = mkpath("plots/L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi]")
    savefig("$folder/Magnetic structure factor L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1] mu = [$current_mu]  phi = [$phi].pdf")

end