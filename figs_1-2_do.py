""" Numerical evaluation of integral in the effective cross section (Haug, AA 191, 181 (1989), eq. (27)) or of the energy loss normalized to density (Haug, eq. 15).  
Here we decompose the integral in the sum of a part with an integrable singularity and a part with a coordinate transformation which eliminates the algebraic part, leaving only an exponential (see paper). This should be numerically more stable. 
We use the mpmath library for high precision quadrature (http://mpmath.org).
We also plot the asynptotic expansion of the integral for tau << 1

THIS VERSION ONLY PRINTS OUTPUTS TO A FILE, AVAILABLE FOR PLOTTING (*_do = data only).

Usage: python sigma_eff_do.py <input parameter file>

Format of the input parameter file:

Electron density (tuple, float, cm^-3)
Temperature (tuple, float, K)
eps_min, eps_max, n_eps (tuple [float, float, integer] - Dimensionless electron total energy (in units of m_e c^2 -> eps > 1), number of eps points for plots of mean free path l_mfp(eps)
maximum degree for integral precision (integer)
"temp" or "ne": Compute either sigma_eff (using "temp" as parameter) or (1/n_e)*de/dt (using "ne" as parameter)
Name of the output data file (string)

"""
import math, sys, ast
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import constants, special
from aux_macros import pause, read_input_pf, gamma_int, z_int, c_eps, \
    sigma_eff, f1, f2, linestyle_str, color
import csv

# Initial setting for mpmath calls.
mp.dps = 35
mp.pretty=True
"""
if(sys.version_info.major < 3):
    plot_par = raw_input('electron density (ne) or temperature (temp) as parameter? ')
else:
    plot_par = input('electron density (ne) or temperature (temp) as parameter? ')
"""

#
# Parameters of the input files
inp_par = read_input_pf(sys.argv[1])

# Use ast.literal_eval to read input from the parameter file into tuple(s).

n_e = ast.literal_eval(inp_par[0].rstrip())  # e^- density in the target (cm^-3)
Tb = ast.literal_eval(inp_par[1].rstrip())    # Target's temperatures (K)
eps_lim = ast.literal_eval(inp_par[2].rstrip()) # Format: eps_min, eps_max, npoints.
max_tsh_degree = int(inp_par[3].rstrip())   # Maximum degree of the tanh-sinh method before leaving
plot_par = inp_par[4].rstrip()
fname_out = inp_par[5].rstrip()     # Output file name

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
log_eps_min = log10(eps_lim[0])
log_eps_max = log10(eps_lim[1])
# Last element has negative index.
n_p = int(eps_lim[-1])

delta_log_eps = (log_eps_max-log_eps_min)/(n_p-1)

eps_i = c_eps(log_eps_min, log_eps_max, delta_log_eps, n_p)

f_out = open(fname_out, 'w')

# Bounds for gamma_r

#gg_1 = 1.1
#gg_2 = 200.

if (plot_par == 'temp'):
    f_out.write("n_rows="+str(n_p)+"\n")
    ne = n_e[0]
    f_out.write("ne (cm^-3)="+str(ne)+"\n")

    
    for temp in Tb:
        f_out.write("T(K) = "+str(temp))
        gamma_r_max = mp.mpf('inf')
        tau, gamma_r_min = gamma_int(temp, ne)
        xy = gamma_r_min
        z_min, z_max = z_int(ne, temp, gamma_r_max)
        sig_mfp_eff = sigma_eff(tau, eps_i, gamma_r_min, z_min, z_max, ne, max_tsh_degree)
        sig_eff = sig_mfp_eff[0]     # Effective cross section (eq. 27 of Haug 1998 - exact)
        mfp = sig_mfp_eff[1]     # Mean free path (eq. 24 and 27 of Haug 1998 - exact)

# Now plots the asymptotic expansion fit.
# We only plot one curve because the fit is independent ot temp.
        sigma_eff_asymp = []
        c_as = (8.*mp.pi)*mp.power(r0_v,2)
        
        for eps in eps_i:
            eps2 = eps*eps
            eps3 = eps*eps*eps
            p2 = (eps2 - 1.)
            p = mp.sqrt(p2)
            aux1 = eps
            aux2 = p*(eps - 1.)
            aux = aux1/aux2
            aux_den1 = tau*(2.*eps3 + 7.*eps2 + 4.*eps + 2.)
            aux_den2 = 2.*p2*eps2
            aux_den = aux_den1/aux_den2
            sigma_eff_asymp_p = c_as*aux*(1. + aux_den)

#
            sigma_eff_asymp.append(sigma_eff_asymp_p)
            
        print(np.array(eps_i,dtype=np.float64).size, np.array(sig_eff,dtype=np.float64).size, np.array(sigma_eff_asymp, dtype=np.float64).size,"\n")
        f_out.write("\n")
        np.savetxt(f_out,np.column_stack((np.array(eps_i,dtype=np.float64), np.array(sig_eff,dtype=np.float64), np.array(sigma_eff_asymp, dtype=np.float64))))
        # ARRIVED HERE (06/03/2024) - Continue by making directly plots.

elif (plot_par == 'ne'):
    f_out.write("n_rows="+str(n_p)+"\n")
    temp = Tb[0]
    f_out.write("Tg (K)="+str(temp)+"\n")

    for ne in n_e:
        f_out.write("ne (cm^-3)="+str(ne))
        f_out.write("\n")
        gamma_r_max = mp.mpf('inf')
        tau, gamma_r_min = gamma_int(temp, ne)
        z_min, z_max = z_int(ne, temp, gamma_r_max)
        dedt = en_loss(tau, eps_i, om_p_pref, \
                       gamma_r_min, z_min, z_max, ne, max_tsh_degree)
        c_as = (2.*mp.pi)*mp.power(r0_v,2)*(cl*1.e2)*ne
        dedt_asymp = []
        for eps in eps_i:
            p = mp.sqrt(eps*eps - 1.)
            aux_num = eps*(tau - (eps - 1.))
            aux_den = p*(eps - 1.)
            dedt_asymp_p = c_as*(aux_num/aux_den)*f_sintheta(eps, om_p_pref, ne)
            dedt_asymp.append(dedt_asymp_p)
        np.savetxt(f_out,np.column_stack((np.array(eps_i,dtype=np.float64), np.array(dedt,dtype=np.float64), np.array(dedt_asymp, dtype=np.float64))))
    
else:
    print("Input should be either 'temp' or 'ne' ")
    sys.exit(0)


