% OP_F_V_BD_TILES: assemble the Neumann contribution to the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) from the trimmed elements.
%
%   rhs = op_f_v_bd_tiles (space, msh_trimmed, iside, coeff);
%
% INPUT:
%
%   space: object representing the space of trial functions (see sp_trimming)
%   msh:   object defining the domain partition and the quadrature rule (see msh_trimming)
%   side: label of the boundary side  where to compute the boundary integral
%   coeff: function handle to compute the Neumann function
%
% OUTPUT:
%
%   rhs:    assembled rhs vector
%
% Copyright (C) 2022 Pablo Antolin, Ondine Chanon, Luca Coradello, Riccardo Puppi
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

function rhs = op_f_v_bd_tiles(space, msh_trimmed, side, f) 

msh_tiles_bd = msh_evaluate_boundary_tiles (msh_trimmed, side, true);
sp_tiles = sp_evaluate_tiles (space, msh_tiles_bd, 'value', true);

for idim = 1:msh_tiles_bd.rdim
  x{idim} = reshape (msh_tiles_bd.geo_map(idim,:,:), msh_tiles_bd.nqn, msh_tiles_bd.nel);
end
coeffs = f (x{:});

rhs = op_f_v (sp_tiles, msh_tiles_bd, coeffs);

end
