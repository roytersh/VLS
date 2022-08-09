% *******************************************************************
% This function computes plasma dielectric tensor
% parameters:  w: complex frequency, kpar: the value of kpar,
% kper: the value of kper, 
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
function y = epsilon(w,kpar,kper,params)

  err = 0;
  nsp = params.nsp;

  y = [1 0 0; 0 1 0; 0 0 1];

  for s=1:nsp
    y = y + chi(w,kpar,kper,s,params,err);
  end
    
end
