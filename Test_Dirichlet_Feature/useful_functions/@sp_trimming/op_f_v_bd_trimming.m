% OP_F_V_BD_TRIMMING: assemble the Neumann contribution to the right-hand side vector r = [r(i)], with  r(i) = (f, v_i) on a boundary of a trimmed domain.
%
%   rhs = op_f_v_bd_trimming (spv, msh_trimmed, side, coeff);
%
% INPUT:
%     
%   spv:       object representing the function space (see sp_trimming)
%   msh:       object defining the domain partition and the quadrature rule (see msh_trimming)
%   side: label of the boundary to compute the boundary integral
%   coeff:     function handle to compute the Neumann contribution
%
% OUTPUT:
%
%   rhs: assembled right-hand side from the Neumann contribution
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

function rhs = op_f_v_bd_trimming(space, msh_trimmed, side, f)

rhs = zeros (space.ndof, 1);

bnd_id = msh_trimming_boundary_find(msh_trimmed.reparam, side);


% Assembly of non trimmed elements
non_trimmed_elems_bd_ids = msh_trimmed.reparam.boundaries(bnd_id).non_trimmed_elem_bd_ids;
nb_non_trimmed_elements = msh_trimmed.reparam.boundaries(bnd_id).nb_non_trimmed_elements;
if nb_non_trimmed_elements > 0
    param_side = msh_trimmed.reparam.boundaries(bnd_id).parametric_side;
    local_indices_current_side_non_trimmed = get_boundary_indices(param_side, msh_trimmed.msh_cart.nel_dir, non_trimmed_elems_bd_ids);
    
    %structure containing the quadrature rule in the untrimmed elements of the physical domain
    msh_elems = msh_eval_boundary_side(msh_trimmed.msh_cart, param_side, local_indices_current_side_non_trimmed); 
    sp = sp_evaluate_element_list(space.space_untrimmed.boundary(param_side), msh_elems);

    x = cell (msh_trimmed.rdim, 1);
    for idim = 1:msh_elems.rdim
        x{idim} = reshape (msh_elems.geo_map(idim,:,:), msh_elems.nqn, msh_elems.nel);
    end
    coeffs = f (x{:});

    dofs = space.space_untrimmed.boundary(param_side).dofs;
    [~, active_dofs, active_dofs_on_boundary] = intersect (space.active_dofs, dofs);
    rhs_loc = op_f_v (sp, msh_elems, coeffs);
    rhs(active_dofs) = rhs(active_dofs) + rhs_loc(active_dofs_on_boundary); % assembling of non-trimmed elements
end


% Assembly of trimmed elements
nb_reparam_elements = msh_trimmed.reparam.boundaries(bnd_id).nb_reparam_elements;
if nb_reparam_elements > 0
    rhs = rhs + op_f_v_bd_tiles(space, msh_trimmed, side, f); 
end

end


