% SP_TRIMMING: constructor of B-spline spaces on trimmed domains.
%
%     sp_trimmed = sp_trimming (knots, degree, msh_trimmed)
%
% INPUTS:
%
%     knots:     open knot vector (cell array of size [1, ndim])
%     degree:    b-spline polynomial degree (vector of size [1, ndim])
%     msh:       msh object that defines the quadrature rule (see msh_trimming)
%   
% OUTPUT:
%
%    sp_trimmed: object of the class representing the discrete function space.
%
%  FIELDS
%  space_untrimmed: the underlying structured B-spline space (see sp_scalar).
%  ncomp: number of components
%  ndof:  number of functions intersecting the trimmed domain (active functions).
%  active_dofs:  indices of all those functions in the original B-spline space.
%  trimmed_dofs: indices of functions that do not vanish on trimmed elements.
%  global_to_active: from global (unstructured) indexing, to active indexing.
%
% Copyright (C) 2021 Luca Coradello, Ondine Chanon, Rafael Vazquez
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

function sp_trimmed = sp_trimming (knots, degree, msh_trimmed)

space = sp_bspline (knots, degree, msh_trimmed.msh_cart);

trimmed_dofs = sp_get_basis_functions (space, msh_trimmed.msh_cart, msh_trimmed.reparam.trimmed_elem_ids);
other_dofs = sp_get_basis_functions (space, msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids);
active_dofs = union (trimmed_dofs, other_dofs);
ndof = numel (active_dofs);
global_to_active = zeros (space.ndof, 1);
global_to_active(active_dofs) = 1:ndof;

sp_trimmed.space_untrimmed = space;
sp_trimmed.ncomp = space.ncomp;
sp_trimmed.ndof = ndof;
sp_trimmed.active_dofs = active_dofs;
sp_trimmed.trimmed_dofs = trimmed_dofs;
sp_trimmed.global_to_active = global_to_active;

sp_trimmed = class (sp_trimmed, 'sp_trimming');

end