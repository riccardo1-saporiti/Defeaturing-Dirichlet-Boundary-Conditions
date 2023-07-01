% SP_GET_BASIS_FUNCTIONS: Compute the indices of trimmed B-splines acting on a list of cells.
%
% [fun_indices, indices_per_cell] = sp_get_basis_functions (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_trimming)
%    msh:     object defining the domain partition and the quadrature rule (see msh_trimming)
%    indices: indices of the cells.
%
% OUTPUT:
%    fun_indices: indices of the basis functions acting on the cells.
%    indices_per_cell: cell-array with indices of the basis functions on each cell.
%
% Copyright (C) 2015, 2016, 2017, 2021 Rafael Vazquez
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

function [function_indices, indices_per_cell] = sp_get_basis_functions (space, msh, cell_indices)

[function_indices, indices_per_cell] = sp_get_basis_functions (space.space_untrimmed, msh.msh_cart, cell_indices);

function_indices = space.global_to_active(function_indices);
if (nargout == 2)
  indices_per_cell = cellfun (@(x) space.global_to_active(x), indices_per_cell, 'UniformOutput', false);
end

end
