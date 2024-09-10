""" Auxiliary functions """
import os, sys, math, string, ast
import numpy as np
import mpmath as mp
from mpmath import *
from scipy import constants
from scipy.special import gamma
#import OpenGL.GL
#import cyglfw3

# This function defined on-the-fly truncates any object containing floats (including arrays) to have 3 decimal digits.
float_formatter = "{:.3f}".format
exp_f = "{:e}".format

# Initial setting for mpmath calls.
mp.dps = 15
mp.pretty=True

linestyle_str = [
     'solid',     # Same as (0, ()) or '-'
     'dotted',    # Same as (0, (1, 1)) or ':'
     'dashed',    # Same as '--'
     'dashdot']  # Same as '-.'

"""
linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or ':'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'
"""
color = ['black',
         'blue',
         'red',
         'green',
         'magenta']

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")
#    programPause = raw_input("Press the <ENTER> key to continue...")

def read_input_pf(inp_par_file):
    usage='Usage: %s Input file' % sys.argv[0]
    out_values = []
    try:
       inp_par_file = sys.argv[1]
    except:
        print(usage); sys.exit[1]
# Opens the input file
    ifile = open(inp_par_file, 'r')
# Reads all the lines. The input file should contain the following lines (in the same order): dir name, prefix name, serial number, serial number,...
# NOTE: no empty lines (neither at the end)!
    out_values = ifile.readlines()
    return out_values

if __name__=='__main__':
    print("Module test_rf works")
#
def k_nu(n, x):
# Approximation to second kind Bessel function of order n (not necessarily integer), from Palade & Pomarjanschi, Romanian Journal of Physics 68, 108 (2023), arxiv:2303.13400.:2303.13400
    pi = mp.pi
    two_nu = 2.*n
    c_nu = (0.2168 + 0.932*n)/(0.392+n)
    aux_1 = mp.power(two_nu,c_nu)
    gam_nu = 2.*aux_1/(1+aux_1)
    aux_2 = mp.gamma(1./gam_nu)
    lam_nu = gam_nu*mp.sqrt(pi)*((mp.gamma(n+0.5)/(mp.gamma(n)*aux_2)))
    arg_knu_1 = mp.power((x/lam_nu),gam_nu)
    k_nu_app = mp.gamma(n)*mp.power(2.,(n-1))/(mp.exp(arg_knu_1)*mp.power(x,n))
    return k_nu_app

#
def gamma_int(Temp, n_e):
# Computes the extremes of the integration interval for gamma.
#
# All quantities in constants are in SI units.
    r0 = constants.physical_constants['classical electron radius']
# tau: K_boltzmann*T_c/(m_e*c^2)
# convert everything to cgs
    k_b = constants.value(u'Boltzmann constant')*1.e7
    m_e = constants.value('electron mass')*1.e3
    cl = constants.value('speed of light in vacuum')*1.e2
    r0_v = constants.value('classical electron radius')*1.e2
    tau_f = k_b/(m_e*(cl**2))
#
    tau = tau_f*Temp

    q = mpf(8.22e-19*n_e/Temp)  # Haug (1989), eq. (29)
    q2 = mp.power(q,2)
    gamma_r_min = mpf(1. + 0.25*q + 0.0625*q2)  # Haug (1989), eq. (28)
    return tau, gamma_r_min

#
def z_int(n_e, Temp, gam_max):
# Interval extremes for z, the transformed variable for the second integral.
# We compute explicitly z_min in terms of q because if the use the expression 
# for z = (gamma - 1/(gamma + 1) the numerator is always =0, due to the
# limited precision of the arithmetic in Python.
#
    q = mpf(8.22e-19*n_e/Temp)  # Haug (1989), eq. (29)
    q2 = mp.power(q,2)
    z_min = (0.25*q + 0.0625*(q**2))/(mpf(2.e0) + 0.25*q + 0.0625*q2)
    if (gam_max == mp.mpf('inf')):
        z_max = 1.
    else:
        z_max = (gam_max - 1.)/(gam_max + 1)

#    z_max = (gam_max - 1.)/(gam_max + 1)
    return z_min, z_max

#
def c_eps(lg_eps_min, lg_eps_max, delta, n_points):
# Fills the tuple of eps1, the dimensionless energy of the projectile particle.
    eps_i=[]
    i = 0
    while i < n_points:
        lg_eps = lg_eps_min + i*delta
        eps_i.append(math.pow(10,lg_eps))
        i += 1
    return eps_i

def sigma_eff(tau, eps_i, g_min, z_min, z_max, n_e, m_deg):
    from scipy import constants, special
# Computes the effective cross section (Haug 1989, eq. 25)
    #
    pi = mp.pi
    r0_v = constants.value('classical electron radius')*1.e2   # (cm)
    cm2pc = mpf(1./3.085678e18)   # cm -> parsec
    sigma_e = []
    l_mean_fp = []
    inv_tau = mp.mpmathify(1./tau)
    
    for eps in eps_i:
#        print("eps=", eps)
        exp_ext = inv_tau*eps
        mp_quad1 = mp.quad(lambda x: f1(x, eps, tau), [g_min, mp.mpf('inf')], maxdegree=m_deg, error=True)
        mp_quad2 = mp.quad(lambda z: f2(z, eps, tau), [z_min, z_max], maxdegree=m_deg, error=True)
        p1 = mp.sqrt(eps*eps - 1)
        sigma_eff_integr = mp.exp(-exp_ext)*(mp_quad1[0]+mp_quad2[0])
#        delta_sigma_eff_integr = mp.exp(-exp_ext)*(mp_quad1[1]+mp_quad2[1])
        sigma_eff = (4.*pi*mp.power(r0_v,2)/(eps*p1))*sigma_eff_integr/mp.besselk(2,inv_tau)
#        delta_sigma_eff = (4.*pi*mp.power(r0_v,2)/(eps*p1))*delta_sigma_eff_integr/mp.besselk(2,inv_tau)
        sigma_e.append(sigma_eff)
#        delta_sigma_e.append(delta_sigma_eff)
        lmfp = mpf((p1/(eps*n_e*sigma_eff))*cm2pc)

        l_mean_fp.append(lmfp)

    return sigma_e, l_mean_fp

#
def en_loss(tau, eps_i, om_p_pf, g_min, z_min, z_max, n_e, m_deg):
# Computes the r.h.s of Haug (1988), eq. (15).
    pi = mp.pi
    r0_v = constants.value('classical electron radius')*1.e2   # (cm)
    cm2pc = mpf(1./3.085678e18)   # cm -> parsec
    en_loss_rhs = []
    inv_tau = mp.mpmathify(1./tau)
    fac_en_l = pi*mp.power(r0_v,2)*(constants.speed_of_light*1.e2)  # light speed here in cgs
    for eps in eps_i:
        p1 = mp.sqrt(eps*eps - 1)
        exp_ext = inv_tau*eps
        mp_quad1 = mp.quad(lambda x: f1_en(x, eps, om_p_pf, n_e, tau), [g_min, mp.mpf('inf')], maxdegree=m_deg, error=True)
        mp_quad2 = mp.quad(lambda z: f2_en(z, eps, om_p_pf, n_e, tau), [z_min, z_max], maxdegree=m_deg, error=True)
        fac_en_l2 = fac_en_l*n_e*mp.exp(-exp_ext)/(eps*p1*mp.besselk(2,inv_tau))
        en_l = fac_en_l2*(mp_quad1[0] + mp_quad2[0])
#        en_l = fac_en_l2*mp_quad2[0]
        en_loss_rhs.append(en_l)
        
    return en_loss_rhs

#
def f1(x, eps1, tau):
#def f1(x, eps1, tau, arg_max):
# The integrand here is exp(+eps_1/tau) times eq. (27) of Haug (1989), thus the argument of the exponential in the integrand becomes: exp(-eps_1*(gamma_r - 1)/tau).
    mp.dps = 15
    mp.pretty=True
    x=mp.mpmathify(x)
    eps1=mp.mpmathify(eps1)
    tau=mp.mpmathify(tau)
    one = mpf(1.e0)
    eps_d = mp.fdiv(eps1,tau)
    eps1_2 = mp.power(eps1,2)
#    c_up = arg_max
#    eps_d = mp.mpmathify(eps1/tau)
#    eps_d = eps1/tau
    eps_d = mp.fdiv(eps1,tau)
    p1 = mp.sqrt(eps1_2 - 1.)
    arg_exp1 = (eps1*(x-1.) - p1*mp.sqrt(x*x - 1.))/tau
    arg_exp2 = (eps1*(x-1.) + p1*mp.sqrt(x*x - 1.))/tau
#    integr_1 = mp.sqrt((x+1)/(x-1))
    integr_1 = mp.sqrt(mp.fdiv((x+one),(x-one)))
#    integr_2 = mp.exp(-arg_exp1) - mp.exp(-arg_exp2)
    integr_2 = mp.fsub(mp.exp(-arg_exp1), mp.exp(-arg_exp2))
#    return integr_1*integr_2
    return mp.fmul(integr_1, integr_2)
#
def f2(z, eps1, tau):
    mp.dps = 15
    mp.pretty = True
    z = mp.mpmathify(z)
    eps1 = mp.mpmathify(eps1)
    tau = mp.mpmathify(tau)
    one = mpf(1.e0)
#
    eps_d = mp.fdiv(eps1,tau)
    x = (1 + z)/(1 - z)
    x2 = mp.power(x,2)
    p1 = mp.sqrt(eps1*eps1 - 1.)
    arg_exp1 = (eps1*(x-1) - p1*mp.sqrt(x2 - 1))/tau
    arg_exp2 = (eps1*(x-1) + p1*mp.sqrt(x2 - 1))/tau
#    integr_1 = 1./(2*mp.power(z,1.5))
    integr_1 = mp.fdiv(one,(2*mp.power(z,1.5)))
#    integr_2 = mp.exp(-arg_g1) - mp.exp(-arg_g2)
    integr_2 = mp.fsub(mp.exp(-arg_exp1), mp.exp(-arg_exp2))
    return mp.fmul(integr_1,integr_2)
#
def f1_en(x, eps1, om_p_pf, ne, tau):
    mp.dps=25
    mp.pretty=True
    x=mp.mpmathify(x)
    eps1=mp.mpmathify(eps1)
    tau=mp.mpmathify(tau)
    one = mpf(1.e0)
    xm1 = mp.fsub(x,one)
    eps_d = mp.fdiv(eps1,tau)
    eps1_2 = mp.power(eps1,2)
    eps_d = mp.fdiv(eps1,tau)
    p1 = mp.sqrt(eps1_2 - 1.)
    x2 = mp.sqrt(x*x - 1.)
    arg_exp1 = (eps1*(x-1.) - p1*x2)/tau
    arg_exp2 = (eps1*(x-1.) + p1*x2)/tau
#
    integr_1a = mp.sqrt(mp.fdiv((x+one),(x-one)))*(eps1*(x - 1.) + tau - p1*x2)*f_sintheta(x, om_p_pf, ne)
    integr_1 = integr_1a*mp.exp(-arg_exp1)
    integr_2a = mp.sqrt(mp.fdiv((x+one),(x-one)))*(eps1*(x - 1.) + tau + p1*x2)*f_sintheta(x, om_p_pf, ne)
    integr_2 = integr_2a*mp.exp(-arg_exp2)
    return mp.fsub(integr_1, integr_2)
#
def f2_en(z, eps1, om_p_pf, ne, tau):
    mp.dps=25
    mp.pretty=True
    z=mp.mpmathify(z)
    eps1=mp.mpmathify(eps1)
    tau=mp.mpmathify(tau)
    one = mpf(1.e0)
    x = (1 + z)/(1 - z)
    eps_d = mp.fdiv(eps1,tau)
    eps1_2 = mp.power(eps1,2)
    eps_d = mp.fdiv(eps1,tau)
    p1 = mp.sqrt(eps1_2 - 1.)
    x2 = mp.sqrt(x*x - 1.)
    arg_exp1 = (eps1*(x-1.) - p1*x2)/tau
    arg_exp2 = (eps1*(x-1.) + p1*x2)/tau
#
    integr_1a = mp.fdiv(one,(2*mp.power(z,1.5)))*(eps1*(x - 1.) + tau - p1*x2)*f_sintheta(x, om_p_pf, ne)
    integr_1 = integr_1a*mp.exp(-arg_exp1)
    integr_2a = mp.fdiv(one,(2*mp.power(z,1.5)))*(eps1*(x - 1.) + tau + p1*x2)*f_sintheta(x, om_p_pf, ne)
    integr_2 = integr_2a*mp.exp(-arg_exp2)
    return mp.fsub(integr_1, integr_2)
#
def theta_mstar(xx, om_p_pf, ne):
    from scipy import constants
    mp.dps = 25
# om_p_pf -> omega_p prefactor, defined in pl_en_loss
    ne *= 1.e6  # Convert density from cm^-3 to m^-3 (SI units)
# om_p: 
    om_p = om_p_pf*mp.sqrt(ne)
    x = mp.mpmathify(xx)
    one = mpf(1.0)
    x2 = mp.fsub(mp.power(x,2),one)
    beta_star = mp.fdiv(mp.sqrt(x2),x, prec=20)
    d = beta_star*constants.speed_of_light/om_p
    if((constants.m_e*constants.speed_of_light*d) <= 0):
        print('test th_mstar', constants.m_e, constants.speed_of_light, x, xx, beta_star, om_p, d, (constants.m_e*constants.speed_of_light*d))
        sys.exit(0)

    theta_1 = 4.*constants.hbar/(constants.m_e*constants.speed_of_light*d)
    theta_2 = x/((x+1)*mp.sqrt(2*(x-1)))
    theta = theta_1*theta_2
    return theta
#
def f_sintheta(xx, om_p_pf, ne):
    import mpmath as mp
    mp.dps = 25
    xx = mp.mpmathify(xx)
    f_sth = 1 - mp.ln(2*mp.power(mp.sin(theta_mstar(xx, om_p_pf, ne)),2))
    f_sth += mp.power(((xx-1.)/xx),2)*(mp.ln(2.) + 1.25e-1)
#    print('inside f_sintheta',theta_mstar(xx, om_p_pf, ne), mp.sin(theta_mstar(xx, om_p_pf, ne)), f_sth, om_p_pf)
    return f_sth
#
def dedt(xx, om_p_pf, c_as, ne, tau):
    import mpmath as mp
    mp.dps = 25
    c_as_ne = c_as*ne     # light speed cl from mks to cgs
    xx = mp.mpmathify(xx)
    p = mp.sqrt(xx*xx - 1.)
    aux_num = xx*(tau - (xx - 1.))
    aux_den = p*(xx - 1.)
    dedt_asymp_p = c_as_ne*(aux_num/aux_den)*f_sintheta(xx, om_p_pf, ne)
    return dedt_asymp_p

#
def er_int(xx, b_j, tau):
# Internal energy of a jet moving with dimensionless speed b_j(=v_jet/c).
# This is the integrand to be used in the numerical evaluation with mpmath.quad.# In units of mc^2
    import mpmath as mp
    b_j2 = mp.power(b_j,2.)
    g_j = 1./mp.sqrt(1. - b_j2)
    lamb = g_j/tau
    a_exp_1 = xx + b_j*sqrt(mp.power(xx,2.)-1)
    arg_exp = -lamb*a_exp_1
    er_in = mp.power(xx,2.)*mp.exp(arg_exp)
    return er_in
#
def ind_match_lines(fname, reg_match):
    ind=[]
    i = 0
    with open(fname, 'r') as f_search:
        for num, line in enumerate(f_search):
            if reg_match in line.split():
                ind.append(num)
                i += 1

    f_search.close()
    return ind

                
