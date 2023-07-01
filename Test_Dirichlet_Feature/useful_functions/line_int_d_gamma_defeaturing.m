% Compute the line integral I = \int_{\gamma} f ds.
%
%   int = line_int_d_gamma_defeaturing (spv, msh, coeff);
%
% INPUT:
%
%   spv:   structure representing the function space (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   coeff: source function evaluated at the quadrature points
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function int = line_int_d_gamma_defeaturing (spv, msh, coeff)

 int = 0;
 
 coeff = reshape (coeff, spv.ncomp, msh.nqn, msh.nel); %rimuovere
 for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
     int = int + sum( msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .* ...
                sum(1 .* coeff(:, :, iel), 1).'); 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'line_int_d_gamma_defeaturing: singular map in element number %d', iel)
    end
 end
end
