% SP_GET_CELLS: Compute the indices of the cells within the support of a list of trimmed B-spline functions.
%
% [cell_indices, indices_per_function] = sp_get_cells (space, msh, indices)
%
% INPUT:
%    space:   object defining the space of discrete functions (see sp_trimming)
%    msh:     object defining the domain partition and the quadrature rule (see msh_trimming)
%    indices: indices of the functions.
%
% OUTPUT:
%    cell_indices: indices of the cells within the support of the basis functions.
%    indices_per_function: indices of the cells within the support of each basis function.
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

function [cell_indices, indices_per_function] = sp_get_cells (space, msh, fun_indices)

fun_indices = space.active_dofs(fun_indices);
[cell_indices, indices_per_function] = sp_get_cells (space.space_untrimmed, msh.msh_cart, fun_indices);

active_cells = union (msh.reparam.trimmed_elem_ids, msh.reparam.non_trimmed_elem_ids);
cell_indices = intersect (cell_indices, active_cells);
indices_per_function = cellfun (@(x) intersect (active_cells, x), indices_per_function, 'UniformOutput', false);

end
