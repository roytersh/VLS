'''
Vlasov linear solver for magnetized plasma

written by Vadim Roytershteyn, SSI, 2025

adapted from an earlier Octave code

'''

import numpy as np
from scipy.special import ive
from scipy.special import wofz

class Vlasov:

    BesselEps = 1E-15
    BesselMaxN = 5000
    root_eps = 1E-10
    maxit = 100
    kpar_min = 1E-8
    sqrt_pi = np.sqrt(np.pi)
    
    def __init__(self,nsp,Tper,Tpar,Upar,ms,Wcs,Wps,verbose):

        # number of species
        self.nsp = nsp
        # T_\perp
        self.Tper = Tper
        # T_||
        self.Tpar = Tpar
        # mass
        self.ms   = ms
        # gyrofreuency [inlcuding sign!]
        self.Wcs  = Wcs
        # plasma frequency
        self.Wps = Wps
        # parallel drift velocity
        self.Upar = Upar  

        # debug flag
        self.verbose = verbose
        
# ----------------------------------------------------------
# This function computes plasma dielectric tensor
# see Stix, Waves in Plasmas, 1992, Chapter 6
# arguments:
# w: complex frequency,
# kpar: the value of kpar,
# kper: the value of kper
# Return value: complex epsilon        
# ----------------------------------------------------------
    def epsilon(self,w,kpar,kper):

        y = np.eye(3,dtype=np.complex128)

        for s in range(self.nsp):
            val,err = self.chi(w,kpar,kper,s)
            y += val
            
        return y

# -------------------------------------------------------------------------    
# This function computes plasma succeptability for maxwellian distributions
# parameters:  w: complex frequency, kpar: the value of kpar,
# kper: the value of kper, s: species index
# see Stix, Waves in Plasmas, 1992, Chapter 6    
# -------------------------------------------------------------------------
    def chi(self,w,kpar,kper,s):

        # relative error for bessel function sum
        eps = self.BesselEps

        # maximum number of Bessel functions
        maxn = self.BesselMaxN

        # T_\perp
        Tper = self.Tper[s]
        # T_||
        Tpar = self.Tpar[s]
        # mass
        ms   = self.ms[s]
        # gyrofreuency [inlcuding sign!]
        Wcs  = self.Wcs[s]
        # plasma frequency
        Wps  = self.Wps[s]
        # parallel drift velocity for species s
        Upar = self.Upar[s]  

        # perp. thermal speed
        vper2 = 2*Tper/ms       
        # perp. thermal speed
        vpar  = np.sqrt(2*Tpar/ms)      

        # argument of bessel functions, etc
        lamb = kper**2*vper2/(2*Wcs**2)  
   
        # contribution from the parallel drift
        y = np.zeros((3,3),dtype=np.complex128)

        if (abs(kpar) > self.kpar_min): 
            y[2,2] += 2*Wps**2*Upar/(w*kpar*vper2)
        else:
            print("Please implement the limit kpar=0")
            exit()
   
        # now add the sum over 
        for nm in range(maxn):
        
           Yn = self.get_yn(nm,w,kpar,kper,lamb,vpar,vper2,Upar,Tper,Tpar,Wcs,Wps)
           dely = Yn*Wps**2/w
           y = y + dely

           err = np.max( np.abs(dely/y) )

           if (err < eps):
               break
    
        # check convergence
        if (err > eps):
            print('Bessel sum did not converge!')


        return y,err

    # ---------------------------
    # plasma dispersion function
    # ---------------------------
    def Z(self, xi):    
        Z_ = 1j*self.sqrt_pi*wofz(xi)
        return Z_

    # ---------------------------------------------
    # This function generates matrices Yn, 
    # see Stix, Waves in Plasmas, 1992, Chapter 6
    # ---------------------------------------------    
    def get_yn(self,nm,w,kpar,kper,lamb,vpar,vper2,Upar,Tper,Tpar,Wcs,Wps):
  
        if (nm==0):
            nmin = 0
        else:
            nmin = -2

        Yn = np.zeros((3,3),dtype=np.complex128)

        # !!!!!!!!!!!!!!!!!!!
        # TESTING ONLY!!!!
        # nmin = 0;
        # !!!!!!!!!!!!!!!!!!!

        for signn in np.arange(1,nmin,-2):
    
            n = nm*signn

            xi = (w - kpar*Upar - n*Wcs)/(kpar*vpar)
            Z0 = self.Z(xi)
            An = 1/w*(Tper-Tpar)/Tpar + Z0*( (w-kpar*Upar -n *Wcs)*Tper +n*Wcs*Tpar)/(kpar*vpar*w*Tpar)
            Bn = ( (w-n *Wcs)*Tper - (kpar*Upar - n*Wcs)*Tpar)/(kpar*w*Tpar) \
                + (w-n *Wcs)/(kpar**2*vpar*w*Tpar)*( (w-kpar*Upar -n *Wcs)*Tper + n*Wcs*Tpar)*Z0

            bopt = 0;
    
            In = ive(n, lamb)      # besseli(n,lamb,1)
            In_p1 = ive(n+1, lamb) # besseli(n+1,lamb,1)
            In_m1 = ive(n-1, lamb) # besseli(n-1,lamb,1)
            dIn_dz = 0.5*(In_m1 + In_p1)
    
            y11 = n**2*In/lamb*An
            y12 = -1j*n*(In-dIn_dz)*An
            y13 = kper/Wcs*n*In*Bn/lamb
            y22 = (n**2/lamb*In+2*lamb*(In-dIn_dz))*An
            y23 = 1j*kper/Wcs*(In-dIn_dz)*Bn
            y33 = 2*(w-n*Wcs)*In*Bn/(kpar*vper2)
       
            Yn += np.array([[y11,y12,y13],[-y12,y22,y23],[y13,-y23,y33]])
       
        return Yn     

    # ---------------------------------------------------------
    # This function returns n x ( n x E) matrix
    # ---------------------------------------------------------
    def M1(self,w,kpar,kper):

        # convert k's to c*k/w
        npar = kpar/w
        nper = kper/w

        # n x ( n x E) matrix
        res = np.array([[-npar**2, 0, npar*nper],[0, -npar**2-nper**2, 0],[npar*nper, 0, -nper**2]])

        return res
    # ---------------------------------------------------------
    # This function computes determinant of n x ( n x E) matrix
    # ---------------------------------------------------------
    def detM(self,w,kpar,kper):


        # add plasma response
        M = self.M1(w,kpar,kper) + self.epsilon(w,kpar,kper)
        
        # compute its determinant
        res = np.linalg.det( M)

        return res

    # ---------------------------------------------------------
    # This function computes eigenvector E & B
    # ---------------------------------------------------------
    def eienvectors(self,w,kpar,kper):

        # first, reduced matrx
        mt = self.M1(w,kpar,kper) + self.epsilon(w,kpar,kper)
        rhs = -mt[:2,0]
        mr  =  mt[:2,1:3]
        er  = np.linalg.solve(mr, rhs)

        # form full vector e and normalize
        esol = np.ones(3,dtype=np.complex128)
        esol[1:3] = er
        
        esol = esol/np.max(abs(esol))

        # find B
        bsol = esol.copy()
        bsol[0] = -kpar*esol[1]/w
        bsol[1] =  (kpar*esol[0] - kper*esol[2])/w
        bsol[2] = kper*esol[1]/w

        return esol,bsol

    # ---------------------------------------------------------
    # rudimentary newton solver
    # arguments: w0(guess), kpar, kper
    # returns w, cflag (convergence)
    # ---------------------------------------------------------
    def my_newton_solve(self,w0,kpar,kper):

        # solver parameters: precision of the root
        root_eps = self.root_eps
        # solver parameters: maximum number of iterations for root finding
        maxit    = self.maxit
        
        it = 0
        w = w0 + 0j
        wstep = 1E-6*w
        err = 1.0
        cflag = True
        
        while (( err > root_eps ) and (it < maxit)):

            # find derivative
            d1 = self.detM(w,kpar,kper)
            d2 = self.detM(w + wstep,kpar,kper)
            dd_dw = (d2-d1)/wstep
            dw = - d1/dd_dw
            w = w + dw
            err = abs(dw/w)
            it += 1

        if (err > root_eps):
            print("Error: Newton iterations did not converge!")
            cflag = False

        return w,cflag

    # ---------------------------------------------------------
    # solve interface
    # ---------------------------------------------------------
    def solve(self,w0,kpar,kper):
        w,cflag = self.my_newton_solve(w0,kpar,kper)
        return w, cflag
