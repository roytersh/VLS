# *******************************************************************
# This code solves linear dispersion relation in bi-Maxwellian plasma
# See e.g. Stix, Waves in Plasmas, 1992, Chapter 6
#
# The code works in units of c,wpe,me,e
#
# Copyright (C) 2025 Vadim Roytershteyn
#
# *******************************************************************

import numpy as np
from vlasov import Vlasov
import matplotlib.pyplot as plt

sqrt = np.sqrt
degree = 1.0/180.0*np.pi

# --------------------------
# specify plasma parameters
# --------------------------

# wpe/wce
wpewce = 4
# mi/me
mime = 1836

# denisities
ne_1 = 0.8        # cold
ne_2 = 1 - ne_1   # hot

# parallel electron beta. NB: these betas are defined with the TOTAL density
beta_ec_par = 6.2500e-04
# perpendicular electron beta
beta_ec_per = beta_ec_par
# parallel ion beta
beta_i_par = beta_ec_par
# perpendicular ion beta
beta_i_per = beta_ec_par

# HOT populations
beta_eh_par = 200*beta_ec_par
# perpendicular electron beta
beta_eh_per = 5*beta_eh_par

# define convenience parameters
wpe = 1.0                   # electron plasma frequency
wpi = 1.0/np.sqrt(mime)     # ion plasma frequency
Wce  = - 1.0/wpewce         # electron cyclotron frequency (with sign!)
Wci  = - Wce/mime           # ion cyclotron frequency

# arrays of plasma parameters for all species
# plasma frequency
Wps = [wpe*sqrt(ne_1), wpe*sqrt(ne_2), wpi]
# cycltron frequency
Wcs = [Wce, Wce, Wci]
# mass ratios
ms  = [1,1,mime]
# parallel drift speed (with offset bulk electron population to make
# sure that the total current is zero 
Upar = [0, 0 , 0]
# parallel in perpendicular temperatures beta=2*n0*T/B0^2
# so T = 0.5*B0^2*beta; in the selected units B0 = 1/wpewce
bf = 0.5/wpewce**2
Tpar = [beta_ec_par*bf, beta_eh_par*bf, beta_i_par*bf]
Tper = [beta_ec_per*bf, beta_eh_per*bf, beta_i_per*bf]

# put all of the parameters into a structure to be passed to sub-routines
nsp = 3

# create solver object
verbose = True
solver = Vlasov(nsp,Tper,Tpar,Upar,ms,Wcs,Wps, verbose)

# ******************************************************************
#  Begin dispersion relation
# ******************************************************************

# in this example, compute w(k) for a given angle theta w.r.t. magnetic field
theta_degree = 1E-3
theta = degree*theta_degree

# values of k to compute w for

nk = 100

# k = np.linspace( 6, 1, nk ) # part A
k = np.linspace( 0.1, 2, nk ) # part B

# pre-compute k_par and k_per
kpar = k*np.cos(theta)
kper = k*np.sin(theta)


E = np.zeros( (nk,3), dtype=np.complex128)
B = np.zeros( (nk,3), dtype=np.complex128)
w = np.zeros( (nk,), dtype=np.complex128)

# ------------------------------------------------------------------
# the guess value for the first k (remember, this is in wpe)
# this is an extremely important parameter that determines
# what branch we analyze . Solver will not converge (or will converge
# to something strange, unless a good guess is given 

wg0 = k[0]**2*np.abs(Wce)

# ------------------------------------------------------------------

# initialize guess value used for each iteration
wg = wg0

# iterate over the wavenumbers
for ik in range(nk):

  kpar_ = kpar[ik]
  kper_ = kper[ik]

  # use newton's iterations to find the root
  w_, cflag = solver.solve(wg, kpar_, kper_)

  # print a status message
  print("ik=%i, k=%g, w=%g, g=%g" % (ik,k[ik],w_.real,w_.imag))

  # next guess  
  wg = w_
  # save the result
  w[ik] = w_

  # find eigenvectors and save the result
  esol,bsol = solver.eienvectors(w_,kpar_,kper_)
  E[ik,:] = esol
  B[ik,:] = bsol


# ********************************************
# save the dispersion relation
# ********************************************

data = np.stack([kpar,kper,w.real,w.imag]).T
for v in range(3):
  tmp = np.stack([ E[:,v].real, E[:,v].imag, B[:,v].real, B[:,v].imag]).T
  print(tmp.shape)
  print(data.shape)  
  data = np.hstack((data,tmp))
fn = 'theta_%g.dat' % theta_degree
hdr = 'kpar, kper, wr, gamma Re(Ex) Im(Ex) Re(Bx) Im(Bx) ... Im(Bz)'
np.savetxt(fn,data,header=hdr)

# ********************************************
# plot the results
# ********************************************
fig, ax = plt.subplots(2,1,sharex=True)

kp = kpar
wp = w.real
gp = w.imag

axc = ax[0]
axc.plot(kp,wp)
axc.set_ylabel(r'$\omega/\omega_{pe}$')

axc = ax[1]
axc.plot(kp,gp)
axc.set_ylabel(r'$\gamma/\omega_{pe}$')

axc.set_xlabel(r'$k{d_e}$')

fn = 'disp.pdf'
fig.savefig(fn,bbox_inches='tight')

plt.show()

#EBratio = sqrt(sum(abs(E(:,:)).^2, 2)./sum(abs(B(:,:)).^2, 2));

