% OP_GRADU_GRADV_TRIMMING: assemble the matrix K = [k(i,j)], k(i,j) = (epsilon grad u_j, grad v_i), for a trimmed geometry.
%
%   mat = op_gradu_gradv_trimming (spu, spv, msh_trimmed, [epsilon]);
%   [rows, cols, values] = op_gradu_gradv_trimming (spu, spv, msh_trimmed, [epsilon]);
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
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
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

function varargout = op_gradu_gradv_trimming(space1, space2, msh_trimmed, coeff)

if (nargin < 4)
  coeff = @(varargin) ones(size(varargin{1}));
end

% Assembly of untrimmed elements
A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);
msh_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids);
sp1 = sp_evaluate_element_list(space1.space_untrimmed, msh_elems, 'gradient', true);
sp2 = sp_evaluate_element_list(space2.space_untrimmed, msh_elems, 'gradient', true);

if(~isempty(msh_trimmed.reparam.non_trimmed_elem_ids))  
    for idim = 1:msh_trimmed.rdim
        x{idim} = reshape (msh_elems.geo_map(idim,:,:), msh_elems.nqn, msh_elems.nel);
    end
    coeffs = coeff (x{:});
    A_loc = op_gradu_gradv (sp1, sp2, msh_elems, coeffs); % assembling of non-trimmed elements
    A = A + A_loc(space2.active_dofs, space1.active_dofs);
end

% Assembly of trimmed elements
if (~isempty(msh_trimmed.reparam.trimmed_elems))
  A = A + op_gradu_gradv_tiles(space1, space2, msh_trimmed, coeff); 
end

if (nargout == 1)
  varargout{1} = A;
elseif (nargout == 3)
  [rows, cols, vals] = find (A);
  varargout{1} = rows;
  varargout{2} = cols;
  varargout{3} = vals;
end

end
