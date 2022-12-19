# VLS
A collection of linear Vlasov solvers in various approximations


Most of the practitioners working in plasma physics end up developing a number of linear solvers. I intend to collect the solvers that I have developed over the years in this repository, in the hope that they will be useful to other people. The codes are described below.


EM_magnetized_octave/

This directory contains a solver for linear dispersion relation in magnetized bi-Maxwellian plasma. The governing equations are presented in many textbooks. The code mostly follows the style and notations of Stix, Waves in Plasmas, 1992, Chapter 6

The code works in units of speed of light c, plasma frequency wpe, electron mass me, and electron mass e. It can be easily adopted for other units. 
The code has been tested with GNU Octave, but should also work in Matlab, perhaps with minor modifications.

The code could be used to analyze equilibria where the distribution function for each plasma species is a Maxwellian (potentially with a drift along the background magnetic field). Arbitrary number of species could be included, which allows one to model a rather large class of distribution functions. Parameters of the distributions are specified in warm_plasma.m 

The code works by finding the root (in the complex plane) of the determinant that arises from linearizing Vlasov-Maxwell system (details could be found in many textbooks, e.g. in the Stix reference above). At present, a very simple Newton root-finding algorithm is used, which requires a rather good initial guess to converge. Note that sometimes it is advantageous to start the search in a region of the dispersion relation where the mode of interest is clearly separated from other neighboring modes and then trace the dispersion relation back to the parameters of interest. Currently, warm_plasma.m is setup to analyze stability of an electron beam propagating nearly parallel to the background. The initial guess is given at relatively large values of the wavenumber and the dispersion is traced back to the intersection between the beam mode and the plasma oscillation, where an instability is found. 
