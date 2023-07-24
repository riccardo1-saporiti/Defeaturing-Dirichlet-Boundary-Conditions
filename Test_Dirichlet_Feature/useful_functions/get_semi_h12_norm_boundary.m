function semi_h12_boundary_norm = get_semi_h12_norm_boundary(problem_data , method_data , alpha )

%Compute the $H^{1/2}$ seminorm of $d_{\gamma}$ over the boundary $\gamma$

problem_data.h = @( x , y , ind ) 1  + alpha * log( sqrt( x .^2 + y.^ 2 ) ) / ( 2*pi)  ; %d_\gamma: difference
%of Dirichlet data "1" and exact solution of the problem with the forcing
%Dirac 

[geometry, msh, space, u] = ...
              mp_solve_laplace (problem_data, method_data);

% vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
% sp_plot_solution ( u , space , geometry , vtk_pts)
problem_data.uex     = @(x, y) 0 * x .* y ;
problem_data.graduex = @(x, y) cat (1, ...
                reshape ( 0 * x .* y , [1, size(x)]), ...
                reshape ( 0 * x .* y , [1, size(x)]));
[~ , ~ , semi_h12_boundary_norm] = sp_h1_error (space, msh, u, problem_data.uex, problem_data.graduex);

semi_h12_boundary_norm = semi_h12_boundary_norm ^ 2 ; 
return