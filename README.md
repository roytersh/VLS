# VLS
A collection of linear Vlasov solvers in various approximations


Most of the practitioners working in plasma physics end up developing a number of linear solvers. I indend to collect the solvers that I have developed over the years in this repository, in the hope that they will be useful to other people.


EM_magnetized_octave/

This directory contains a solver for linear dispersion relation in bi-Maxwellian plasma. The governing equations are presented in many textbooks. The code mostly follows the style and notations of Stix, Waves in Plasmas, 1992, Chapter 6

The code works in units of spped of light c, plasma frequency wpe, electron mass me, and electron mass e. It can be easily adopted for other units. 
The code has been tested with GNU Octave, but should also work in Matlab, perhaps with minor modifications.
