include("utils.jl")

# Load data.
ground_states, L, t, U_0, U_1, mu, phi = load_data(
    "output/(0.1) L = [12]  t = [1]  U_0 = [0.1]  U_1 = [0.05]  phi = [0.5235987755982988].h5")


# Real space lattice vectors.
a = [1, 0]
b = [1/2, 1/2]

# Reciprocal lattice vectors.
#R = [0 -1; 1 0] # 90 degree rotation matrix.
#a_reciprocal = (2*pi*R*b)/(dot(a, R*b))
#b_reciprocal = (2*pi*R*a)/(dot(b, R*a))

site_unit_cell_coords = [[0, 0]]

N_x = 12
N_y = 2

x, y = lattice_coords_2d(a, b, site_unit_cell_coords, N_x, N_y)
# x_reciprocal, y_reciprocal = lattice_coords_2d(a_reciprocal, b_reciprocal, site_unit_cell_coords, N_x, N_y)

######################### --PLOTS-- ####################################

plot_density_vs_mu(ground_states, L, t, U_0, U_1, mu, phi)

#plot_density_vs_mu_all_U("output", L, t, phi)

mu_index = 21

plot_lattice_occupation(ground_states, L, t, U_0, U_1, mu, phi, mu_index, x, y)

plot_spin_correlations(ground_states, L, t, U_0, U_1, mu, phi, mu_index, x, y)

# Charge structure factor - 1d FT of Ntot.
plot_charge_structure_factor(ground_states, L, t, U_0, U_1, mu, phi, mu_index)

# Magnetic structure factor - 1d FT of spin_corr.
plot_magnetic_structure_factor(ground_states, L, t, U_0, U_1, mu, phi, mu_index)