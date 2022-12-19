% *******************************************************************
% This code solves linear dispersion relation in bi-Maxwellian plasma
% See e.g. Stix, Waves in Plasmas, 1992, Chapter 6
%
% The code works in units of c,wpe,me,e
% It has been tested with GNU Octave,
% but should also work in Matlab, perhaps with minor modifications
%
% Copyright (C) 2022 Vadim Roytershteyn
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
% *******************************************************************

% clear the workspace
clear all;

% solver parameters: precision of the root
root_eps = 1E-8;
% solver parameters: maximum number of iterations for root finding
maxit    = 30;

% --------------------------
% specify plasma parameters
% --------------------------

% wpe/wce
wpewce = 1.75;
% mi/me
mime = 4*1836;
% parallel electron beta. NB: these betas are defined with the TOTAL density
beta_e_par = 2.46e-6;
% perpendicular electron beta
beta_e_per = beta_e_par;
% parallel ion beta
beta_i_par = beta_e_par;
% perpendicular ion beta
beta_i_per = beta_e_par;

% in this exampe, we add an electron beam with speed Vb 
Vb = 0.27;
% the beam has a guassian distribution with temperature Tb
Tb = (0.1*Vb)**2;
% beam density
nb = 1E-2;

% define convenience parameters
wpe = 1.0;                   % electron plasma frequency
wpi = 1/sqrt(mime);          % ion plasma frequency
Wce  = - 1/wpewce;           % electron cyclotron frequency (with sign!)
Wci  = - Wce/mime;           % ion cyclotron frequency

% arrays of plasma parameters for all species
% plasma frequency
Wps = [wpe wpi*sqrt(1+nb) wpe*sqrt(nb)];
% cycltron frequency
Wcs = [Wce Wci Wce];
% mass ratios
ms  = [1  mime 1];
% parallel drift speed (with offset bulk electron population to make
% sure that the total current is zero 
Upar = [ -nb*Vb 0 Vb];
% parallel in perpendicular temperatures beta=2*n0*T/B0^2
% so T = 0.5*B0^2*beta; in the selected units B0 = 1/wpewce
bf = 0.5/wpewce^2;
Tpar = [beta_e_par*bf beta_i_par*bf Tb];
Tper = [beta_e_per*bf beta_i_per*bf Tb];

% put all of the parameters into a structure to be passed to sub-routines
params.nsp = 3;
params.Tpar = Tpar;
params.Tper = Tper;
params.Wps = Wps;
params.ms = ms;
params.Wc = Wcs;
params.Upar = Upar;


% ******************************************************************
%  Begin dispersion relation
% ******************************************************************

% in this example, compute w(k) for a given angle theta w.r.t. magnetic field
theta = 1/180.0*pi;
% values of k to compute w for
k = linspace( 10, 3, 100 );
% pre-compute k_par and k_per
kpar = k*cos(theta);
kper = k*sin(theta);
nk = length(kper);

% ------------------------------------------------------------------
% the guess value for the first k (remember, this is in wpe)
% this is an extremely important parameter that determines
% what branch we analyze . Solver will not converge (or will converge
% to something strange, unless a good guess is given 
wg0 = 1.0001
  % ------------------------------------------------------------------

% initialize guess value used for each iteration
wg = wg0;

% iterate over the wavenumbers
for ik = 1:nk

  kpar_ = kpar(ik);
  kper_ = kper(ik);

  % convert k's to c*k/w
  npar = @(w) kpar_./w;
  nper = @(w) kper_./w;

  % n x ( n x E) matrix
  M1 = @(w) [-npar(w).^2 0 npar(w).*nper(w); 0 -npar(w).^2-nper(w).^2 0;  npar(w).*nper(w) 0 -nper(w).^2];

  % form full matrix and compute its determinant
  M = @(w) det( epsilon(w,kpar_,kper_,params) + M1(w));

  % use newton's iterations to find the root

  it = 0;
  w = wg;
  wstep = 1E-6*w;

  do

     % find derivative
     d1 = M(w);
     d2 = M(w + wstep);
     dd_dw = (d2-d1)/wstep;
     dw = - d1/dd_dw;
     w = w + dw;
     err = abs(dw/w);
     it++;
  until ((err < root_eps)||(it == maxit))

  if (err > root_eps)
    printf("Error: Newton iterations did not converge!\n");
  endif

  % print a status message
  printf("ik=%i, k=%g, w=%g, g=%g\n",ik,k(ik),real(w),imag(w)); 
  fflush(stdout);

  % next guess  
  wg = w;
  % save the result
  ws(ik) = w;

 % compute eigenvector

  % first, reduced matrx
  mt = epsilon(w,kpar_,kper_,params) + M1(w);
  rhs = -mt(1:2,1);
  mr  =  mt(1:2,2:3);
  er  = mr\rhs;

  % form full vector e and normalize
  esol = [1; er];
  esol = esol/max(abs(esol));

  % find B
  bsol(1) = -kpar_*esol(2)/w;
  bsol(2) =  (kpar_*esol(1) - kper_*esol(3))/w;
  bsol(3) = kper_*esol(2)/w;

  % save the result
  E(ik,:) = esol;
  B(ik,:) = bsol;
  
endfor

% ********************************************
% save the dispersion relation
% ********************************************

data = [k.' real(ws.') imag(ws.') ];
save -ascii "dispersion.dat" data

% ********************************************
% plot the results
% ********************************************
figure(1)
clf; hold on;
plot(k,real(ws),'-r','linewidth',5)
xlabel('k_{d_e}')
ylabel('\omega/\omega_{pe}')
print -dpdf -F:18 "omega.pdf"

figure(2)
clf; hold on;
plot(k,imag(ws),'-r','linewidth',5)
xlabel('k{d_e}')
ylabel('\gamma/\omega_{pe}')
print -dpdf -F:18 "gamma.pdf"

EBratio = sqrt(sum(abs(E(:,:)).^2, 2)./sum(abs(B(:,:)).^2, 2));
figure(3)
clf; hold on;
semilogy(k,EBratio,'-r','linewidth',5)
xlabel('k{d_e}')
ylabel('|E|/|B|')
print -dpdf -F:18 "comp.pdf"

