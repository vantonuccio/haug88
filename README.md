# haug88
We make here available Python macros and auxiliary files connected to our paper on the calculation of Haug (1988) cross sections and electron loss rate.

Python version >= 3.6 required.

The "Asympototic expansion" pdf file describes with some detail both the asymptotic expansion and the numerical integration scheme adopted to evaluate the coros section ad the energy loss rate.

Reference: Haug, E., AA 191, 181 (1988) (Haug 1988).

## REQUIREMENTS:

Libraries: numpy, scipy.
We make use of the mpmath library (https://mpmath.org/) for arbitrary precision mathematics, more specifically to compute Haug 1988 exact quadrature formulas for both quantities (cross section and energy loss). 

### List of functions:<br>
**figs_1-2_do.py**: Computes either sigma_eff (using "temp" as parameter) or (1/n_e)*de/dt (using "ne" as parameter), and outputs to an ascii file with (among others) three columns: (energy, exact value, symp. exp.), where "exact value" refers to the numerical quadrature evaluation of the exact expressions from Haug 1988.<br>
*Input*: An ascii file containing parameters.<br>
*Output*: If using "temp" as parameter: ASCII file, three columns, n_eps lines (eps, sigma_eff_exact, sigma_eff_asymp) for each input value of temperature. If using "ne" as parameter: ASCII file, two columns, n_eps lines (eps, de/dt) for each input value of target's electron density.

### List of classes:<br>

haug88.py: Container of a series of functions computing electrons cross sections and energy losses.<br>

Here the list of functions:<br>
**sigma_eff** : Evaluates cross section using eq. (25) of Haug 1988. "quad" Numerical quadrature function from mpmath is used.<br>
*Input*: eps_i, a list of energy values.<br>
*Output*: sigma_eff (list), mean free path (list).<br>

**f1, f2**: The integrands used in sigma_eff.<br>

**sig_eff_as**: Asymptotic expansion of eq. (25) of Haug 1988.<br>
*Input*: eps_i, a list of energy values. <br>
*Output*: sigma_eff (list).<br>

**en_loss**: Computes the r.h.s. of Haug (1988), eq. (15).<br>
*Input*: eps_i, a list of energy values.<br>
*Output*: energy loss (list).<br>

**f1_en, f2_en**: The integrands used in en_loss.<br>
**theta_mstar**: eq. 21a of Haug 1988.<br>
**f_sintheta**: eq 12 of Haug 1988.<br>

**dedt**: Asymptotic expansion of Haug (1988), eq. 15.<br>
*Input*: eps_i, a list of energy values.<br>
*Output*: energy loss (list).
