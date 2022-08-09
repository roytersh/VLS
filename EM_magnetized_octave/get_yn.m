% *******************************************************************
% This function generates matrices Yn, 
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

function Yn = get_yn(nm,w,kpar,kper,lambda,vpar,vper2,Upar,Tper,Tpar,Wcs,Wps)
  
  if (nm==0)
    nmin = 1;
  else
    nmin = -1;
  endif

  Yn = zeros(3,3);

  # !!!!!!!!!!!!!!!!!!!
  # TESTING ONLY!!!!
  # nmin = 1;
  # !!!!!!!!!!!!!!!!!!!

  for signn = 1:-2:nmin
    
    n = nm*signn;

    xi = (w - kpar*Upar - n*Wcs)/(kpar*vpar);
    Z0 = 1i*sqrt(pi)*exp(-xi^2) - 2*dawson(xi);
    An = 1./w*(Tper-Tpar)./Tpar + Z0*( (w-kpar*Upar -n *Wcs)*Tper +n*Wcs*Tpar)./(kpar*vpar*w*Tpar);
    Bn = ( (w-n *Wcs)*Tper - (kpar*Upar - n*Wcs)*Tpar)./(kpar*w*Tpar) ...
	 + (w-n *Wcs)./(kpar^2*vpar*w*Tpar)*( (w-kpar*Upar -n *Wcs)*Tper + n*Wcs*Tpar)*Z0;

    bopt = 0;
    
    In = besseli(n,lambda,1);
    In_p1 = besseli(n+1,lambda,1);
    In_m1 = besseli(n-1,lambda,1);
    dIn_dz = 0.5*(In_m1 + In_p1);
    
    y11 = n^2*In/lambda*An;
    y12 = -1i*n*(In-dIn_dz)*An;
    y13 = kper/Wcs*n*In*Bn/lambda;
    y22 = (n^2/lambda*In+2*lambda*(In-dIn_dz))*An;
    y23 = 1i*kper/Wcs*(In-dIn_dz)*Bn;
    y33 = 2*(w-n*Wcs)*In*Bn/(kpar*vper2);
       
    Yn = Yn + [y11 y12 y13; -y12 y22 y23; y13 -y23 y33];
       
  endfor

end
