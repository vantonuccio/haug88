# haug88
We make here available Python macros and auxiliary files connected to our paper on the calculation of Haug (1988) cross sections and electron loss rate.

Python version >= 3.6 required.

The "Asympototic expansion" pdf file describes with some detail both the asymptotic expansion and the numerical integration scheme adopted to evaluate the coros section ad the energy loss rate.

REQUIREMENTS:

Libraries: numpy, scipy.
We make use of the mpmath library (https://mpmath.org/) for arbitrary precision mathematics, more specifically to compute Haug 1988 exact quadrature formulas for both quantities (cross section and energy loss). 

List of functions:

figs_1-2_do.py: Computes either sigma_eff (using "temp" as parameter) or (1/n_e)*de/dt (using "ne" as parameter), and outputs to an ascii file with (among others) three columns: (energy, exact value, symp. exp.), where "exact value" refers to the numerical quadrature evaluation of the exact expressions from Haug 1988.
Input: An ascii file containing parameters.
