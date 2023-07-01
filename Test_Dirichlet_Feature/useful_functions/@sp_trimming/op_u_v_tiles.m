% OP_U_V_TILES: assemble the contribution to the matrix M = [m(i,j)], m(i,j) = (epsilon u_j, v_i) from the trimmed elements.
%
%   mat = op_u_v_tiles (spu, spv, msh_trimmed, [epsilon]);
%
% INPUT:
%
%   spu:     object representing the space of trial functions (see sp_trimming)
%   spv:     object representing the space of test functions (see sp_trimming)
%   msh:     object defining the domain partition and the quadrature rule (see msh_trimming)
%   epsilon: function handle to compute the diffusion coefficient (optional)
%
% OUTPUT:
%
%   mat:    assembled mass matrix
%
% Copyright (C) 2020 Ondine Chanon, Luca Coradello, Riccardo Puppi
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

function M = op_u_v_tiles (space1, space2, msh_trimmed, coeff)

msh_tiles = msh_evaluate_tiles (msh_trimmed);
sp_tiles1 = sp_evaluate_tiles (space1, msh_tiles, 'value', true);
sp_tiles2 = sp_evaluate_tiles (space2, msh_tiles, 'value', true);

if (nargin == 4)
  for idim = 1:msh_trimmed.rdim
    x{idim} = reshape (msh_tiles.geo_map(idim,:,:), msh_tiles.nqn, msh_tiles.nel);
  end
  coeffs = coeff (x{:});
else
  coeffs = ones (current_msh_tile.nqn, current_msh_tile.nel);
end

M = op_u_v (sp_tiles1, sp_tiles2, msh_tiles, coeffs);

end
