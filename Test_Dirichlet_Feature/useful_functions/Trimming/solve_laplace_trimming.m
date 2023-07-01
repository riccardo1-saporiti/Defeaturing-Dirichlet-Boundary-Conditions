% SOLVE_LAPLACE_TRIMMING: Solve a Laplace problem with a B-spline discretization on a trimmed domain. 
%
% The function solves the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega trimmed
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% USAGE:
%
%  [geometry, msh_cart, space, u, stiff_mat, Nitsche_mat] = solve_laplace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:  labels of the boundariesNeumann boundary condition (may be empty)
%    - drchlt_sides: labels of the boundarieswith strong Dirichlet boundary condition
%    - weak_drchlt_sides: (optional) labels of the boundaries with weak Dirichlet boundary
%                         condition, imposed with Nitsche's method.
%    - c_diff:       diffusion coefficient (epsilon in the equation)
%    - f:            source term
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%    - h:            function for Dirichlet boundary condition
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - reparam:    a structure containing the data related to trimming: 
%       - nb_trimmed_elems: number of trimmed elements
%       - nb_non_trimmed_elems: number of non trimmed elements
%       - non_trimmed_elem_ids: IDs of all the elements that are NOT trimmed
%       - trimmed_elem_ids: IDs of all the elements that are trimmed
%       - trimmed_elems: a structure containing for each trimmed elements the info
%                         coming from the reparametrization tool, such as element ID and tiles structure  
%       - boundaries: a structure containing for each boundary the info
%                      related to the trimmed and non trimmed elements of the
%                      associated boundary
%    - Nitsche_type: (if weak_drchlt_sides is not empty) choice of Nitsche's 
%                    method, it can be either 'symmetric' or 'non symmetric'
%    - parameter: (if weak_drchlt_sides is not empty) penalization
%                 parameter for Nitsche's method
%    - stabilization: boolean flag, (if weak_drchlt_sides is not empty)
%       if true the stabilization of [Buffa,Puppi,Vazquez 2019]
%       
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh_cart: cartesian mesh object that defines the quadrature rule (see msh_trimming)
%  space:    space object that defines the discrete space (see sp_trimming)
%  u:        the computed degrees of freedom
%  stiff_mat: the diagonally-scaled resulting stiffness matrix 
%  Nitsche_mat: (if weak_drchlt_sides is not empty) stiffness matrix
%               corrected by Nitsche's penalization term
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

function [geometry, msh_trimmed, sp_trimmed, u, stiff_mat, varargout] = ...
    solve_laplace_trimming (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end
if ~isfield(problem_data, 'weak_drchlt_sides')
    weak_drchlt_sides = [];
end

% Construct geometry structure
geometry  = geo_load(geo_name);
[knots, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);

% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (zeta, rule);
% msh_cart = msh_cartesian (zeta, qn, qw, geometry);
msh_trimmed = msh_trimming (zeta, qn, qw, geometry, reparam);

% Construct space structure
% space = sp_bspline (knots, degree, msh_cart);
sp_trimmed = sp_trimming (knots, degree, msh_trimmed);

% Assemble the matrices
stiff_mat = op_gradu_gradv_trimming (sp_trimmed, sp_trimmed, msh_trimmed);
rhs = op_f_v_trimming (sp_trimmed, msh_trimmed, f);

% Apply Neumann boundary conditions
for side = nmnn_sides
    gside = @(varargin) g(varargin{:}, side);
    rhs = rhs + op_f_v_bd_trimming(sp_trimmed, msh_trimmed, side, gside);
end

% impose weak Dirichlet b.c. via Nitsche
if ~isempty(weak_drchlt_sides)
    if (~exist ('Nitsche_type', 'var'))
      Nitsche_type = 'symmetric';
    end
    if (~exist('stabilization_type', 'var') || ~stabilization)
        [A_Nitsche, rhs_Nitsche, N] = sp_weak_drchlt_bc_laplace_trimming...
            (sp_trimmed, msh_trimmed, weak_drchlt_sides, hfun, Cpen, Nitsche_type);
    else 
        if (~exist ('stabilization_type', 'var') && msh_trimmed.rdim == 2)
            stabilization_type = 'physical';
        elseif (~exist ('stabilization_type', 'var') && msh_trimmed.rdim == 3)
            stabilization_type = 'parametric';
        end
        [A_Nitsche, rhs_Nitsche, N] = sp_weak_drchlt_bc_laplace_trimming_stab(sp_trimmed, msh_trimmed,...
            weak_drchlt_sides, hfun, Cpen, Nitsche_type, stabilization_type, theta);
    end
    
    Nitsche_mat =  stiff_mat + N;
    stiff_mat = stiff_mat - A_Nitsche;
    rhs = rhs + rhs_Nitsche;
end

% Apply strong Dirichlet boundary conditions
for side = drchlt_sides
    bnd_id = msh_trimming_boundary_find(msh_trimmed.reparam, side);
    param_side = msh_trimmed.reparam.boundaries(bnd_id).parametric_side;
    if param_side == 0
        error('Strong Dirichlet boundary conditions must be applied to non trimmed boundaries!');
    end
end

[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp_trimmed, msh_trimmed, hfun, drchlt_sides);
int_dofs = setdiff (1:sp_trimmed.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system

% preconditioner left-right (PAPP^-1x=Pb), always applied for trimmed problem
K_red = stiff_mat(int_dofs, int_dofs);

if isempty(weak_drchlt_sides)
    D = 1./sqrt(diag(K_red));
else
    Nitsche_mat = Nitsche_mat(int_dofs, int_dofs);
    D = 1./sqrt(diag(Nitsche_mat));
    varargout{1} = Nitsche_mat;
end
n = length(D);
P_pc = spdiags(D(:),0,n,n);
C_pc = P_pc*K_red*P_pc;

b = P_pc*rhs(int_dofs);
u = zeros (sp_trimmed.ndof, 1);
u(drchlt_dofs) = u_drchlt;
u(int_dofs) = P_pc*(C_pc\b);

% cond_numb_pc = condest(C_pc);
% fprintf('Condition number after Jacobi preconditioning: %d. \n', cond_numb_pc);
% fprintf('Residual: %d. \n', norm(stiff_mat(int_dofs, int_dofs)*u(int_dofs) - rhs(int_dofs)));
% 
% stiff_mat = K_red;
% 
end

