% OP_F_V_TILES: assemble the contribution to the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) from the trimmed elements.
%
%   rhs = op_f_v_tiles (space, msh_trimmed, coeff);
%
% INPUT:
%
%   space: object representing the space of trial functions (see sp_trimming)
%   msh:   object defining the domain partition and the quadrature rule (see msh_trimming)
%   coeff: function handle to compute the source function
%
% OUTPUT:
%
%   rhs:    assembled rhs vector
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

function rhs = op_f_v_tiles (space, msh_trimmed, coeff)

msh_tiles = msh_evaluate_tiles (msh_trimmed);
sp_tiles = sp_evaluate_tiles (space, msh_tiles, 'value', true);

for idim = 1:msh_trimmed.rdim
  x{idim} = reshape (msh_tiles.geo_map(idim,:,:), msh_tiles.nqn, msh_tiles.nel);
end
coeffs = coeff (x{:});

rhs = op_f_v (sp_tiles, msh_tiles, coeffs);

end