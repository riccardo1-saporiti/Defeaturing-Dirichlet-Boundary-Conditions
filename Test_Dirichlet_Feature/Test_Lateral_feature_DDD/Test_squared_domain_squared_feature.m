clc
clearvars
close all
epsilon_ref = 1e-2;  
eps_values = epsilon_ref./ 2.^(0:9);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 0 * x .* y ; 
theta = @( x , y ) atan2( y , x );
problem_data.h_1 = @(x, y, ind) 0 * x .* y ;
X = [ 0 0.125 0.25 0.75 0.875 1 ];
V = [ 0 0.5 -0.5 -0.5 0.5 1 ];
[ lagrange_interp , grad_function_part_x ] =  lagrange_sym(X,V);
[problem_data, problem_data_0] = determineBC(problem_data);
problem_data.h = @(x, y, ind)  lagrange_interp( x ) .* ( 1 + sin( 2 * pi * ( y - 1 ) / 4 ) ) .* ( isempty( find( ind == problem_data.drchlt_sides_inhomogeneous ) ) == 0  );
problem_data_0.h = @(x, y, ind) lagrange_interp( x ) .* ( 1 + sin( 2 * pi * ( y - 1 ) / 4 ) ) .* ( isempty( find( problem_data_0.drchlt_sides_inhomogeneous == ind ) ) == 0 ); %%cos( theta( x - 0.5 , y - 1 ) );
problem_data_0.h_auxiliar = @(x, y, ind) lagrange_interp( x ) .* ( 1 + sin( 2 * pi * ( y - 1 ) / 4 ) ) .* ( isempty( find( [problem_data_0.gamma_sides{:}] == ind ) ) == 0 ); %%cos( theta( x - 0.5 , y - 1 ) );


perturbation_type = 'regular'; 

ext = 4;
plotIt = 'true' ; 





method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [ 128 128 ]; % [32 32]; 
method_data.nquad = [5 5];


nn = numel(eps_values);
errh1s = zeros(1,nn);
errh1s_interface = zeros(1,nn);
est = zeros(1,nn);
normu = zeros(1,nn);
errh1s_rel = zeros(1,nn);
meas_gamma = zeros(1,nn);


for ii = 1:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
   
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
    
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_defeaturing_ (problem_data, method_data);
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
    figure
    subplot( 1 , 2 , 1 )
    sp_plot_solution (u, space, omega, vtk_pts)

    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_defeaturing_ (problem_data_0, method_data);
%     strm = struct( msh_0 );( 1 - epsilon )
%     for ii = 1:7
%         vect_plot = struct( strm.msh_patch{ii} ).breaks ;
%         sp_plot_solution (u_0 * 0, space_0, omega_0, vect_plot)
%         
%     end
    subplot( 1 , 2 , 2 )
    sp_plot_solution (u_0 , space_0 , omega_0 , vtk_pts)


 
    % 4) COMPUTE ERROR AND ESTIMATOR
    errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
  
    [est( ii ) , est_h1_tg( ii ) ] = estimator_negative_dirichlet(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data_0.h_auxiliar,...
                        problem_data_0.omega_patches, space , epsilon , perturbation_type );
    
end


close all 

%%
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r' , 'LineWidth', 2 );
    hold on 
    loglog(  eps_values, est_h1_tg , '+-b' , 'LineWidth', 2)
    loglog( eps_values , est , '+-y' , 'Linewidth' , 2 )
    % loglog( eps_values , (1 ./ sqrt( abs( log( eps_values ) ) )) * errh1s( 1 ) * sqrt( abs( log( eps_values( 1 )))) , 'Linewidth' , 3)
    loglog( eps_values, eps_values , 'g:' , 'LineWidth' , 5 )
    loglog( eps_values, eps_values .^ (1.8) * errh1s( 1 ) / ( eps_values( 1 ) ^ 1.8 ), 'k:' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$',  '$Estimator(\nabla_{t}u_0)$' , '$Estimator(u_0)$', '$\epsilon$' , '$\epsilon^{1.8}$',  'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$g_{d}:P_{5}(x)(1+sin(\frac{2\pi (y-1)}{4}))\:' ...
        '\partial{\Omega}\setminus{\Gamma}\:\:Lateral\:squared\:section\:square$'] , 'Interpreter','latex' , 'FontSize', 30 ...
    )
%     if saveIt
%         saveas(fig, filename, 'epsc');
%     end
end

ord_eps_est_h1 = log2( est_h1_tg( 1 : end - 1 ) ./ est_h1_tg( 2 : end ) )
ord_eps_est = log2( est( 1 : end - 1 ) ./ est( 2 : end ) )
ord_eps_err = log2( errh1s( 1 : end - 1 ) ./ errh1s( 2 : end ) )










%%
function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    srf_0(1) = nrb4surf ([0 1-epsilon], [1-epsilon 1-epsilon], [0 1], [1-epsilon 1]);
    srf_0(1) = nrbkntins(srf_0(1), {1-ext*epsilon/(2*(1-epsilon)), 0.5 });
    
    srf_0(2) = nrb4surf ([1-epsilon 1-epsilon], [1 1-epsilon], [1-epsilon 1], [1 1]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, 0.5});
    
    srf_0(3) = nrb4surf ([0 0], [1-epsilon 0], [0 1-epsilon], [1-epsilon 1-epsilon]);
    srf_0(3) = nrbkntins(srf_0(3), {1-ext*epsilon/(2*(1-epsilon)), 1-ext*epsilon/(2*(1-epsilon))});
    
    srf_0(4) = nrb4surf ([1-epsilon 0], [1 0], [1-epsilon 1-epsilon], [1 1-epsilon]);
    srf_0(4) = nrbkntins(srf_0(4), {0.5, 1-ext*epsilon/(2*(1-epsilon))});
    
   for ii = 1 : 4
       nrbplot( srf_0(ii ) , [ 20 20 ] )
       hold on 
       view( 2 )
   end

    srf = srf_0([1, 3:4]);
    srf_F = srf_0(2);

end




function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [ ];
    problem_data.drchlt_sides = [ 1 4 5 7 2 8 3 6 ] ;
    problem_data.drchlt_sides_homogeneous = [1 4 5 7 ];
    problem_data.drchlt_sides_inhomogeneous = [2 8 3 6];
    problem_data.gamma_sides = [2 8];
    problem_data.omega0_patches = [ 1:3];
%     problem_data.gamma_info = {[2 3 4], {4, 4, 4}, [1 3; 2 3; 3 3]};

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = [ 1 5 6 8 3 4 2 7 ] ; 
    problem_data_0.drchlt_sides_homogeneous =  [ 1 5 6 8 ] ;
    problem_data_0.drchlt_sides_inhomogeneous = [3 4 2 7];
    problem_data_0.omega_patches = [ 1 3 4 ];
    problem_data_0.gamma_sides = cell(4,1); 
    problem_data_0.gamma_sides(1) = {2}; % relative to each patch
    problem_data_0.gamma_sides(4) = {4}; % relative to each patch
end