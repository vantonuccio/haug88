import math as m
import mpmath as mp
from mpmath import *
import numpy as np
from scipy import constants
from scipy.integrate import quad
from aux_haug88 import z_int


import sys
#print(h5py.__version__)
#print(np.__version__)

# ------------------------------------------------------------------
# Class containing asymptotic expansions to espressions for e^- cross
# sections and energy loss, from Haug A&A 191, 181 (1998).
# 
# ------------------------------------------------------------------

class haug(dict):
    '''
    '''
    def __init__(self, npoints, emin, emax):
        
        '''
        '''
        #
        # Initial setting for mpmath calls.
        mp.dps = 15
        mp.pretty=True
        #
        pi = np.pi

#        import numpy as np
        self.npoints = npoints
        self.emin = emin
        self.emax = emax

        '''
        self.zeta_min = None
        self.zeta = None
        '''
        self.ne = 1.e10
        self.temp = 1.e6
        self.gam_max = 1.e2
        self.zminmax = z_int(self.ne, self.temp, self.gam_max)

        #

    def sigma_eff(self, tau, eps_i, g_min, z_min, z_max, n_e, m_deg):
        '''
Computes the cross section sigma_eff and the mean free path lmfp, eqs. (25) and (24) of Haug, AA 191, 181 (1988) (herefater Haug88), respectively, for a list of input electron's energies eps_i and an input target's electron density n_e.

Uses the quad function of the mpmath library for numerical integration.

Output: lists sigma_e, l_mean_fp.
        '''

        from scipy import constants, special
        # 
        # Few constants 
        pi = mp.pi
        r0_v = constants.value('classical electron radius')*1.e2   # (cm)
        cm2pc = mpf(1./3.085678e18)   # cm -> parsec
        sigma_e = []
        l_mean_fp = []
        inv_tau = mp.mpmathify(1./tau)
        for eps in eps_i:
            #
            exp_ext = inv_tau*eps
            mp_quad1 = mp.quad(lambda x: self.f1(x, eps, tau), [g_min, mp.mpf('inf')], maxdegree=m_deg, error=True)
            mp_quad2 = mp.quad(lambda z: self.f2(z, eps, tau), [z_min, z_max], maxdegree=m_deg, error=True)
            p1 = mp.sqrt(eps*eps - 1)
            sigma_eff_integr = mp.exp(-exp_ext)*(mp_quad1[0]+mp_quad2[0])
            #
            sigma_eff = (4.*pi*mp.power(r0_v,2)/(eps*p1))*sigma_eff_integr/mp.besselk(2,inv_tau)
            #
            sigma_e.append(sigma_eff)
            lmfp = mpf((p1/(eps*n_e*sigma_eff))*cm2pc)
            l_mean_fp.append(lmfp)

        return sigma_e, l_mean_fp

    def f1(self,x, eps1, tau):
# See auxiliary file 
# The integrand here is exp(+eps_1/tau) times eq. (27) of Haug (1989), thus the argument of the exponential in the integrand becomes: exp(-eps_1*(gamma_r - 1)/tau).
        mp.dps = 15
        mp.pretty=True
        x=mp.mpmathify(x)
        eps1=mp.mpmathify(eps1)
        tau=mp.mpmathify(tau)
        one = mpf(1.e0)
        eps_d = mp.fdiv(eps1,tau)
        eps1_2 = mp.power(eps1,2)
        eps_d = mp.fdiv(eps1,tau)
        p1 = mp.sqrt(eps1_2 - 1.)
        arg_exp1 = (eps1*(x-1.) - p1*mp.sqrt(x*x - 1.))/tau
        arg_exp2 = (eps1*(x-1.) + p1*mp.sqrt(x*x - 1.))/tau
        integr_1 = mp.sqrt(mp.fdiv((x+one),(x-one)))
        integr_2 = mp.fsub(mp.exp(-arg_exp1), mp.exp(-arg_exp2))
        return mp.fmul(integr_1, integr_2)

    def f2(self,z, eps1, tau):
        mp.dps = 15
        mp.pretty = True
        z = mp.mpmathify(z)
        eps1 = mp.mpmathify(eps1)
        tau = mp.mpmathify(tau)
        one = mpf(1.e0)
        eps_d = mp.fdiv(eps1,tau)
        x = (1 + z)/(1 - z)
        x2 = mp.power(x,2)
        p1 = mp.sqrt(eps1*eps1 - 1.)
        arg_exp1 = (eps1*(x-1) - p1*mp.sqrt(x2 - 1))/tau
        arg_exp2 = (eps1*(x-1) + p1*mp.sqrt(x2 - 1))/tau
        integr_1 = mp.fdiv(one,(2*mp.power(z,1.5)))
        integr_2 = mp.fsub(mp.exp(-arg_exp1), mp.exp(-arg_exp2))
        return mp.fmul(integr_1,integr_2)

    def sig_eff_as(self,eps_i, tau, c):
        '''
        Asymptotic expansion for sigma_eff (tau << 1).
        We have verified that the for loop is slightly better performing than a series of vectorized numpy array reduce operations. 
        '''

        sigma_eff_asymp = []
        '''
#Numpy vectorized form, slightly worse performing than the for loop
        eps = np.array(eps_i, dtype=np.float64)
        n_eps = eps.size
        eps2 = np.square(eps)
        eps3 = np.power(eps,3.)
        p2 = np.subtract(eps2,1.)
        p = np.sqrt(p2)
        aux1 = eps
        aux2 = np.multiply(p,np.subtract(eps,1.))
        aux = np.divide(aux1,aux2)
        aux_d23 = np.multiply(eps3,2.)
        aux_d72 = np.multiply(eps2,7.)
        aux_d41 = np.multiply(eps,4.)
        aux_dd = [aux_d23, aux_d72, aux_d41]
        aux_d = np.sum(aux_dd)
        aux_den1 = np.multiply(aux_d, tau)
        aux_den2 = np.multiply(2.,np.multiply(p2,eps2))
        aux_den = np.divide(aux_den1,aux_den2)
        sigm_1 = np.multiply(aux, c)
        one = np.ones(n_eps)
        print("test ", one, one.size, aux_den.size)
        sigm_2 = np.add(aux_den, one)
        sigma_eff_asymp_p = np.multiply(sigm_1, sigm_2)
        sigma_eff_asymp.append(sigma_eff_asymp_p)
 
        '''

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
            sigma_eff_asymp_p = c*aux*(1. + aux_den)
            sigma_eff_asymp.append(sigma_eff_asymp_p)

        return sigma_eff_asymp

    def en_loss(self,tau, eps_i, om_p_pf, g_min, z_min, z_max, n_e, m_deg):
        '''
        Computes the electron's energy loss, eq. 15 of Haug 1988.
        '''
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
            en_loss_rhs.append(en_l)

        return en_loss_rhs

    def f1_en(self,x, eps1, om_p_pf, ne, tau):
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

    def f2_en(self, z, eps1, om_p_pf, ne, tau):
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

    def theta_mstar(self, xx, om_p_pf, ne):
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

    def f_sintheta(self, xx, om_p_pf, ne):
        import mpmath as mp
        mp.dps = 25
        xx = mp.mpmathify(xx)
        f_sth = 1 - mp.ln(2*mp.power(mp.sin(theta_mstar(xx, om_p_pf, ne)),2))
        f_sth += mp.power(((xx-1.)/xx),2)*(mp.ln(2.) + 1.25e-1)
#    print('inside f_sintheta',theta_mstar(xx, om_p_pf, ne), mp.sin(theta_mstar(xx, om_p_pf, ne)), f_sth, om_p_pf)
        return f_sth

    def dedt(self, xx, om_p_pf, c_as, ne, tau):
        '''
        Computes the asymptotic expansion of the electron's energy loss.
        '''
        import mpmath as mp
        mp.dps = 25
        c_as_ne = c_as*ne     # light speed cl from mks to cgs
        xx = mp.mpmathify(xx)
        p = mp.sqrt(xx*xx - 1.)
        aux_num = xx*(tau - (xx - 1.))
        aux_den = p*(xx - 1.)
        dedt_asymp_p = c_as_ne*(aux_num/aux_den)*f_sintheta(xx, om_p_pf, ne)
        return dedt_asymp_p

    def er_int(self, xx, b_j, tau):
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
