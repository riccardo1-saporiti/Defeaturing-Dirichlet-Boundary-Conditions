% OP_F_V_TRIMMING: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) in a trimmed domain.
%
%   rhs = op_f_v_trimming (spv, msh_trimmed, coeff);
%
% INPUT:
%     
%   spv:   object representing the function space (see sp_trimming)
%   msh:   object defining the domain partition and the quadrature rule (see msh_trimming)
%   coeff: function handle to compute the source function (optional)
%
% OUTPUT:
%
%   rhs: assembled right-hand side
%
% Copyright (C) 2020 Ondine Chanon, Luca Coradello, Riccardo Puppi
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function rhs = op_f_v_trimming(space, msh_trimmed, coeff)

rhs = zeros (space.ndof, 1);

% Assembly of untrimmed elements
if(~isempty(msh_trimmed.reparam.non_trimmed_elem_ids))
    msh_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids); % structure containing the quadrature rule in the untrimmed elements of the physical domain
    sp = sp_evaluate_element_list(space.space_untrimmed, msh_elems, 'gradient', true); % compute the basis functions in a given list of elements

    for idim = 1:msh_elems.rdim
        x{idim} = reshape (msh_elems.geo_map(idim,:,:), msh_elems.nqn, msh_elems.nel);
    end
    coeffs = coeff (x{:});

    rhs_loc = op_f_v (sp, msh_elems, coeffs);
    rhs = rhs + rhs_loc(space.active_dofs); % assembling of non-trimmed elements
end

% Assembly of trimmed elements
if(~isempty(msh_trimmed.reparam.trimmed_elems))
    rhs = rhs + op_f_v_tiles(space, msh_trimmed, coeff);
end

end
