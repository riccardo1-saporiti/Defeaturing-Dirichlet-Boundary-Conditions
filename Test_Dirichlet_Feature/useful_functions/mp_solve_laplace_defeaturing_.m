 % MP_SOLVE_LAPLACE: solve the Laplacian problem in a multipatch geometry.
%
% Example to solve the diffusion problem
%
%    - div ( epsilon(x) grad (u)) = f    in Omega
%                epsilon(x) du/dn = g    on Gamma_N
%                               u = h    on Gamma_D
%
% where the domain \Omega is formed by several patches of the form F((0,1)^n).
%
% USAGE:
%
%  [geometry, msh, space, u] = 
%          mp_solve_laplace (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
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
%
% OUTPUT:
%
%  geometry: array of geometry structures (see geo_load)
%  msh:      multipatch mesh, consisting of several Cartesian meshes (see msh_multipatch)
%  space:    multipatch space, formed by several tensor product spaces plus the connectivity (see sp_multipatch)
%  u:        the computed degrees of freedom
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2010, 2011, 2013, 2015, 2017 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [geometry, msh, space, u ] = ...
              mp_solve_laplace_defeaturing_ (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% Construct geometry structure, and information for interfaces and boundaries
[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (geo_name);
npatch = numel (geometry);

msh = cell (1, npatch); 
sp = cell (1, npatch);
for iptc = 1:npatch

          [knots{iptc}, zeta{iptc}] = ...
             kntrefine (geometry(iptc).nurbs.knots, nsub-1, degree, regularity);
    
    % Compute the quadrature rule
      rule      = msh_gauss_nodes (nquad);
      [qn, qw]  = msh_set_quad_nodes (zeta{iptc}, rule);
      msh{iptc} = msh_cartesian (zeta{iptc}, qn, qw, geometry(iptc));
    
    % Evaluate the discrete space basis functions in the quadrature points
      sp{iptc} = sp_bspline (knots{iptc}, degree, msh{iptc});
  
end

msh = msh_multipatch (msh, boundaries);

strm = struct( msh );


space = sp_multipatch (sp, msh, interfaces, boundary_interfaces);
clear sp

% Compute and assemble the matrices 
stiff_mat = op_gradu_gradv_mp (space, space, msh, c_diff);
rhs = op_f_v_mp (space, msh, f);

% Add eventual additional Delta term
if ( exist( 'Dirac_Test' , 'var' ) )
    iptc = 1 ;
    %loop over all the patches to find the one that contains the center of
    %the Dirac Delta (could be done one time)
    while  (iptc <= npatch)
    
        [ param_point , conv_indicator ] = nrbinverse( geometry( iptc ).nurbs , dirac_center );
        if conv_indicator
            break 
        end
        iptc = iptc + 1 ;
    
    end
    iptc_star = iptc ;
    qn_delta  = cell( 1 , 2 ) ; 
    qn_delta{ 1 } = param_point( 1 ) ; 
    qn_delta{ 2 } = param_point( 2 ) ; 
    mknots = length (knots{ iptc_star }{ 1 })-1;
    mcp    = -p - 1 + mknots;
    
    s_1 = findspan( mcp , p , qn_delta{ 1 } , knots{iptc_star}{  1}  )  + 1 ;  
    s_2 = findspan( mcp , p , qn_delta{ 2 } , knots{iptc_star}{  2}  )  + 1 ;  
    brk = cell( 1 , 2 ); 
    brk{ 1 } = knots{iptc_star}{ 1 }( s_1 : s_1 + 1 ) ; 
    brk{ 2 } = knots{iptc_star}{ 2 }( s_2 : s_2 + 1 ) ;
    msh_1_pt = msh_cartesian( brk , qn_delta , [] , geometry( iptc_star ) ) ; 
    sp_1_pt = space.sp_patch{ iptc_star}.constructor( msh_1_pt  );  
    msh_prec = msh_precompute( msh_1_pt ); 
    
    if strcmp( singularity_type , 'cos' )
        sp_prec = sp_precompute( sp_1_pt , msh_1_pt , 'value' , true , 'gradient' , true );
        grad_x_shape_fcns = reshape( sp_prec.shape_function_gradients( 1 , : , : ) , [ length( qn_delta{ 1 } ) , sp_1_pt.nsh_max ] );
        rhs_delta = op_f_v_mp_dirac (space, msh, f , iptc_star , - grad_x_shape_fcns , sp_prec.connectivity , alpha_dirac );
        
    elseif strcmp( singularity_type , 'sin' )
        sp_prec = sp_precompute( sp_1_pt , msh_1_pt , 'value' , true , 'gradient' , true );
        grad_y_shape_fcns = reshape( sp_prec.shape_function_gradients( 1 , : , : ) , [ length( qn_delta{ 1 } ) , sp_1_pt.nsh_max ] );
        rhs_delta = op_f_v_mp_dirac (space, msh, f , iptc_star , - grad_y_shape_fcns , sp_prec.connectivity , alpha_dirac );
    elseif strcmp( singularity_type , 'cos_sin' )
        sp_prec = sp_precompute( sp_1_pt , msh_1_pt , 'value' , true , 'gradient' , true );
        grad_x_shape_fcns = reshape( sp_prec.shape_function_gradients( 1 , : , : ) , [ length( qn_delta{ 1 } ) , sp_1_pt.nsh_max ] );
        grad_y_shape_fcns = reshape( sp_prec.shape_function_gradients( 1 , : , : ) , [ length( qn_delta{ 1 } ) , sp_1_pt.nsh_max ] );
        rhs_delta = op_f_v_mp_dirac (space, msh, f , iptc_star , - grad_y_shape_fcns - grad_x_shape_fcns , sp_prec.connectivity , alpha_dirac );
    
    else
        sp_prec = sp_precompute( sp_1_pt , msh_1_pt );
        rhs_delta = op_f_v_mp_dirac (space, msh, f , iptc_star , sp_prec.shape_functions , sp_prec.connectivity , alpha_dirac );
    end
    rhs = rhs + rhs_delta ; 
end

% Apply Neumann boundary conditions
Nbnd = cumsum ([0, boundaries.nsides]);
for iref = nmnn_sides
  iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
  gref = @(varargin) g(varargin{:},iref);
  rhs_nmnn = op_f_v_mp (space.boundary, msh.boundary, gref, iref_patch_list);
  rhs(space.boundary.dofs) = rhs(space.boundary.dofs) + rhs_nmnn;
end

u = zeros (space.ndof, 1);
% Apply Dirichlet boundary conditions 
[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space, msh, h, drchlt_sides);
u(drchlt_dofs) = u_drchlt;

int_dofs = setdiff (1:space.ndof, drchlt_dofs);
rhs(int_dofs) = rhs(int_dofs) - stiff_mat(int_dofs, drchlt_dofs)*u_drchlt;

% Solve the linear system
u(int_dofs) = stiff_mat(int_dofs, int_dofs) \ rhs(int_dofs);

end

%!demo
%! ex_laplace_Lshaped_mp

%!demo
%! ex_laplace_cube_mp

%!demo
%! ex_laplace_thick_L_mp
