% OP_WEAK_DRCHLT_BC_LAPLACE_TRIMMING: assembles the terms needed to impose
%   weak Dirichlet boundary conditions using Nitsche's method on a trimmed
%   geometry.
%
%   [A, rhs, N] = op_weak_drchlt_bc_laplace_trimming (space, msh,...
%       weak_drchlt_sides, hfun, Cpen, Nitsche_type)
%
% INPUT:
%
%   space: object representing the space of trial and test functions (see sp_trimming)
%   msh: object defining the domain partition and the quadrature rule (see msh_trimming)
%   weak_drchlt_sides: labels of boundaries with weak Dirichlet boundary condition
%   hfun: function for Dirichlet boundary condition
%   Cpen: penalization parameter for Nitsche's method
%   Nitsche_type: choice of Nitsche's method, it can be either 'symmetric' or 'non symmetric'
%
% OUTPUT:
%
%   A:    Nitsche correction to the stiffness matrix, dependent on
%         the penalization parameter and on Nitsche_type
%   rhs:  term to be added to the right hand side, coming from Nietsche's method
%   N:    Nitsche penalization matrix N = Cpen*1/h*[(u_j, v_i)]
%
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

function [A, rhs, N] = sp_weak_drchlt_bc_laplace_trimming (space, msh_trimmed,...
    weak_drchlt_sides, hfun, Cpen, Nitsche_type)

if (nargin < 7 || strcmpi(Nitsche_type, 'symmetric'))
    Nitsche_sign = 1;
elseif (strcmpi(Nitsche_type, 'non symmetric'))
    Nitsche_sign = -1;
else
    error ('you must specify if you want symmetric or non symmetric for Nitsche')
end

A = spalloc(space.ndof, space.ndof, 3*space.ndof);
N = spalloc(space.ndof, space.ndof, 3*space.ndof);
rhs = zeros(space.ndof, 1);

if msh_trimmed.reparam.nb_non_trimmed_elems > 0
    mesh_h = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids(1));
    h_el = 1/Cpen*mesh_h.element_size(1);
else
    mesh_h = msh_evaluate_element_list (msh_trimmed.msh_cart, 1);
    h_el = 1/Cpen*mesh_h.element_size(1);
end

act_dofs = space.active_dofs;

for side = weak_drchlt_sides

    bnd_id = msh_trimming_boundary_find (msh_trimmed.reparam, side);

    % Non trimmed part
    non_trimmed_elems_bd_ids = msh_trimmed.reparam.boundaries(bnd_id).non_trimmed_elem_bd_ids;
    if not(isempty(non_trimmed_elems_bd_ids))
        param_side = msh_trimmed.reparam.boundaries(bnd_id).parametric_side;
        local_indices_non_trimmed = get_boundary_indices(param_side, msh_trimmed.msh_cart.nel_dir, non_trimmed_elems_bd_ids);
        mesh_side = msh_eval_boundary_side(msh_trimmed.msh_cart, param_side, local_indices_non_trimmed);
        mesh_side_from_interior = msh_boundary_side_from_interior(msh_trimmed.msh_cart, param_side);
        sp_bnd = space.space_untrimmed.constructor (mesh_side_from_interior);
        mesh_side_from_interior = msh_evaluate_element_list...
            (mesh_side_from_interior, local_indices_non_trimmed);
        sp_bnd = sp_evaluate_element_list (sp_bnd, mesh_side_from_interior,...
            'value', true, 'gradient', true);
        coeff = ones(mesh_side.nqn, mesh_side.nel);

        B = op_gradv_n_u(sp_bnd, sp_bnd, mesh_side, coeff);
        C = op_u_v(sp_bnd, sp_bnd, mesh_side, coeff);
        gradv_n_g = op_gradv_n_f (sp_bnd, mesh_side, ...
            hfun(mesh_side.geo_map(1,:), mesh_side.geo_map(2,:), side));
        g_cdot_v = op_f_v (sp_bnd, mesh_side,...
            hfun(mesh_side.geo_map(1,:), mesh_side.geo_map(2,:), side));

        N = N + 1/h_el*C(act_dofs,act_dofs);
        A = A + B(act_dofs,act_dofs).' + B(act_dofs,act_dofs)*Nitsche_sign - 1/h_el*C(act_dofs,act_dofs);
        rhs = rhs - gradv_n_g(act_dofs)*Nitsche_sign + 1/h_el*g_cdot_v(act_dofs);
    end

    % Reparameterized part
    msh_tiles_bd = msh_evaluate_boundary_tiles (msh_trimmed, side, true);
    sp_tiles = sp_evaluate_tiles (space, msh_tiles_bd, 'value', true, 'gradient', true);

    coeff = ones (msh_tiles_bd.nqn, msh_tiles_bd.nel);

    B = op_gradv_n_u(sp_tiles, sp_tiles, msh_tiles_bd, coeff);
    gradv_n_g = op_gradv_n_f (sp_tiles, msh_tiles_bd,...
        hfun(msh_tiles_bd.geo_map(1,:), msh_tiles_bd.geo_map(2,:), side));
    C = op_u_v(sp_tiles, sp_tiles, msh_tiles_bd, coeff);
    g_cdot_v = op_f_v (sp_tiles, msh_tiles_bd, ...
        hfun(msh_tiles_bd.geo_map(1,:), msh_tiles_bd.geo_map(2,:), side));

    N = N + 1/h_el*C;
    A = A + B.' + B*Nitsche_sign - 1/h_el*C;
    rhs = rhs - gradv_n_g*Nitsche_sign + 1/h_el*g_cdot_v;


end


end
