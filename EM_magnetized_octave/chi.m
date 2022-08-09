% *******************************************************************
% This function computes plasma succeptability for maxwellian distributions
% parameters:  w: complex frequency, kpar: the value of kpar,
% kper: the value of kper, s: species index
% params: structure with the parameters, err: error flag
% see Stix, Waves in Plasmas, 1992, Chapter 6
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
 function y = chi(w,kpar,kper,s,params,errflag)

   eps  = 1E-15; % relative error of Bessel sum calculation
   maxn = 5000; % maximum number of terms in the series

   Tper = params.Tper(s);   % T_\perp
   Tpar = params.Tpar(s);  % T_||
   ms   = params.ms(s);     % mass
   Wcs  = params.Wc(s);     % gyrofreuency (inlcuding sign!)
   Wps  = params.Wps(s);    % plasma frequency
   Upar = params.Upar(s);   % parallel drift velocity for species s
   
   vper2 = 2*Tper/ms;       % perp. thermal speed
   vpar  = sqrt(2*Tpar/ms);      % perp. thermal speed
   
   lambda = kper^2*vper2/(2*Wcs^2);  % argument of bessel functions, etc
   
				% contribution from the parallel drift
   y = [0 0 0; 0 0 0; 0 0 2*Wps^2*Upar./(w*kpar*vper2)];
   
				% now add the sum over n     
   
   
				% error control
   for nm=0:maxn

     Yn = get_yn(nm,w,kpar,kper,lambda,vpar,vper2,Upar,Tper,Tpar,Wcs,Wps);

     y = y + Yn*Wps^2/w  ;

     err = max( abs(Yn(:)./y(:)) );
     if (err < eps) 
       break
     endif
    
   endfor

  if (err > eps)

    printf('Bessel sum did not converge!\n');
    
  endif
  % end sum
       
 end
