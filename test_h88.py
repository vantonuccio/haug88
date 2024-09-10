import math, sys, time, ast  #Note: library ast only available for Python 3, v. >= 3.6
import numpy as np
print("My numpy version is: ", np.__version__)
import mpmath as mp
from mpmath import *
from scipy import constants, special
from aux_haug88 import pause, read_input_pf, gamma_int, z_int, c_eps, \
    sigma_eff, f1, f2, linestyle_str, color
import haug88 as haug88

#
# Parameters of the input files
inp_par = read_input_pf(sys.argv[1])

# Use ast.literal_eval to read input from the parameter file into tuple(s).

n_e = ast.literal_eval(inp_par[0].rstrip())  # e^- density in the target (cm^-3)
Tb = ast.literal_eval(inp_par[1].rstrip())    # Target's temperatures (K)
epsr = ast.literal_eval(inp_par[2].rstrip()) # Format: eps_min, eps_max, npoints
max_tsh_degree = int(inp_par[3].rstrip())   # Maximum degree of the tanh-sinh method before leaving
plot_par = inp_par[4].rstrip()
fname_out = inp_par[5].rstrip()     # Output file name

emin = epsr[0]
emax = epsr[1]
npoints = epsr[-1]

hi = haug88.haug(npoints,emin,emax)


r0_v = constants.value('classical electron radius')*1.e2   # In cgs
r0_v = mpf(r0_v)

cm2pc = mpf(1./3.085678e18)   # cm -> parsec

# This is the constant 4*h_planck_bar/m_e*c from eq (21a) of Haug (1988)
c1 = 4*constants.hbar
e_el = constants.elementary_charge
cl = constants.speed_of_light
c_theta_mstar = c1/(constants.m_e*constants.speed_of_light)

# Prefactor of plasma frequency omega_p (Haug 1988, p. 183, after eq. 21)
om_p_pref = mp.sqrt(4*mp.pi/constants.m_e)*e_el

#
log_eps_min = log10(emin)
log_eps_max = log10(emax)
delta_log_eps = (log_eps_max-log_eps_min)/(npoints-1)

# Creates a list of input energy values, emin <= eps[i] <= emax . The function c_eps assumes equispaced values in logarithmic space.
eps_i = c_eps(log_eps_min, log_eps_max, delta_log_eps, npoints)

f_out = open(fname_out, 'w')

ne = n_e[0]
f_out.write("ne (cm^-3)="+str(ne)+"\n")

# Cross section for temperatures T read from input file, cgs units.

for temp in Tb:
# Monitors execution time.
#    start_time = time.perf_counter()
    f_out.write("T(K) = "+str(temp))
    gamma_r_max = mp.mpf('inf')
    tau, gamma_r_min = gamma_int(temp, ne)
    xy = gamma_r_min
    z_min, z_max = z_int(ne, temp, gamma_r_max)
    sig_mfp_eff = hi.sigma_eff(tau, eps_i, gamma_r_min, z_min, z_max, ne, max_tsh_degree)
    sig_eff = sig_mfp_eff[0]     # Effective cross section (eq. 27 of Haug 1998 - exact)
    mfp = sig_mfp_eff[1]     # Mean free path (eq. 24 and 27 of Haug 1998 - exact)

# Now computes the asymptotic expansion fit of sigma_eff (cgs units)
# We only plotcompute one curve because the fit is independent ot temp.
#
    c_as = (8.*mp.pi)*mp.power(r0_v,2)

    sigma_eff_asymp = hi.sig_eff_as(eps_i, tau, c_as)

#    finish_time = time.perf_counter()
#    print("Program finished in {} seconds".format(finish_time-start_time))

    
    print(np.array(eps_i,dtype=np.float64).size, np.array(sig_eff,dtype=np.float64).size, np.array(sigma_eff_asymp, dtype=np.float64).size,"\n")
    f_out.write("\n")
    np.savetxt(f_out,np.column_stack((np.array(eps_i,dtype=np.float64), np.array(sig_eff,dtype=np.float64), np.array(sigma_eff_asymp, dtype=np.float64))))


