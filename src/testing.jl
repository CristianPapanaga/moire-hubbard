include("utils.jl")
include("moire_hubbard.jl")

###########################################################
#                   DMRG PARAMETERS                       #
# https://itensor.org/docs.cgi?page=tutorials/dmrg_params #
########################################################### 

nsweeps = 3
maxdim = [20, 50, 100]
cutoff = 1E-5

#nsweeps = 5
#maxdim = [20, 50, 100, 150]
#cutoff = 1E-6

###########################################################
#                   SYSTEM PARAMETERS                     #
###########################################################

L = 12
t = 1
mu = -2:0.05:-1
U = 0.1
scale = 0.5
U_0 = U
U_1 = round(U*scale, digits = 5)
phi = pi/6

###########################################################
#                       RUN DMRG                          #
###########################################################

energies, ground_states = moire_hubbard(L, nsweeps, maxdim, cutoff, t, U_0, U_1, mu, phi)

###########################################################
#                       SAVE DATA                         #
###########################################################

# Save groundstates.
# Place any values you later wish to retrieve in [].
# The values in () are used to sort them in natural order, as HDF5 returns datasets in alphabetical order.
let
    file = h5open("output/($U) L = [$L]  t = [$t]  U_0 = [$U_0]  U_1 = [$U_1]  phi = [$phi].h5", "w")
    for (i, psi0) in pairs(ground_states)

        write(file, "($i) mu = [$(mu[i])]", psi0)

    end
    close(file)
end