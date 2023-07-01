clc
clearvars
close all
addpath('useful_functions')
epsilon_ref = 1e-2 ;  
eps_values = epsilon_ref./ 2.^(0:6);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y) 0 * x .* y ; 
theta = @( x , y ) atan2( y , x );

[problem_data, problem_data_0] = determineBC(problem_data);
%impose unitary Dirichlet boundary condition over \gamma and homogeneous
%dirichlet data over \Gamma
problem_data.h = @(x, y, ind) exp( cos( theta( x,  y ) ) ) * ( isempty( find( ind == problem_data.drchlt_sides_inhomogeneous ) ) == 0  ) + 0 * x .* y  ;
problem_data_0.h = @(x, y, ind) 0 * x .* y ;

%auxiliary function used to compute the estimator: to be set equal to the
%dirichlet data over \gamma and non vanishing only on these boundary
problem_data_0.h_auxiliar = @(x, y, ind) exp( cos( theta( x,  y ) ) ) * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;
perturbation_type = 'singular' ; 

ext = 4;
plotIt = 'true' ; 

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [ 8 8 ]; % [32 32]; 
method_data.nquad = [5 5];

nn = numel(eps_values);
est = zeros(1,nn);
est_h1_tg = est ; 
errh1s = est ; 



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
    %%%Code to plot the mesh 
    %     strm = struct( msh_0 );
    %     for ii = 1:7
    %         vect_plot = struct( strm.msh_patch{ii} ).breaks ;
    %         sp_plot_solution (u_0 * 0, space_0, omega_0, vect_plot)
    %         
    %     end
    subplot( 1 , 2 , 2 )
    sp_plot_solution (u_0 , space_0 , omega_0 , vtk_pts)

    % 4) COMPUTE ERROR AND ESTIMATOR
    errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);
    [est( ii ), est_h1_tg( ii ) ] = estimator_negative_dirichlet(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data_0.h_auxiliar ,...
                        problem_data_0.omega_patches, space , epsilon , perturbation_type  );
end


close all 

%%
if plotIt
    fig = figure;
    loglog(eps_values, errh1s , '+-r' , 'LineWidth', 4);
    hold on 
    loglog(  eps_values, est, 'b' , 'Marker', '*' , 'MarkerSize', 10 )
    loglog( eps_values , est_h1_tg , '+-' , 'LineWidth', 2  )
    % loglog( eps_values, eps_values, 'k:' , 'LineWidth' , 5 )
    % loglog( eps_values , 1 ./ sqrt( abs( log(eps_values) ) ) , '+-' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$', 'Estimator', '$Estimator(\nabla_{t}u_0)$' , 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$g_{d}=e^{cos\theta}\:on\:\gamma\:\:\:\:f=0\:in\:\Omega\:\:$'] , 'Interpreter','latex' , 'FontSize', 40 ...
    )

end


ord_eps_est_h1 = log2( est_h1_tg( 1 : end - 1 ) ./ est_h1_tg( 2 : end ) )
ord_eps_est = log2( est( 1 : end - 1 ) ./ est( 2 : end ) )
ord_eps_err = log2( errh1s( 1 : end - 1 ) ./ errh1s( 2 : end ) )
%%

function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    L = sqrt(2)*epsilon/4;
    srf_0(1) = nrbruled (nrbcirc(epsilon, [0 0], pi/4, 3*pi/4), nrbcirc(1, [0 0], pi/4, 3*pi/4));
    srf_0(2) = nrbtform(srf_0(1), vecrotz(pi/2));
    srf_0(3) = nrbtform(srf_0(2), vecrotz(pi/2));
    srf_0(4) = nrbtform(srf_0(3), vecrotz(pi/2));
    for ii = 1:4
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (ext-1)*epsilon/(1-epsilon)});
    end
    srf_0(5) = nrbruled (nrbline([L L], [-L L]), nrbcirc(epsilon, [0 0], pi/4, 3*pi/4));
    srf_0(6) = nrbtform(srf_0(5), vecrotz(pi/2));
    srf_0(7) = nrbtform(srf_0(6), vecrotz(pi/2));
    srf_0(8) = nrbtform(srf_0(7), vecrotz(pi/2));
    for ii = 5:8
        srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L -L], [L, -L], [-L L], [L L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:4);
    srf_F = srf_0(5:9);
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [];
    problem_data.drchlt_sides_inhomogeneous = [1 3 5 7];
    problem_data.drchlt_sides_homogeneous = [2 4 6 8];
    problem_data.drchlt_sides = [ 1 3 5 7 2 4 6 8 ] ;
    problem_data.gamma_sides = [1 3 5 7];
    problem_data.omega0_patches = 1:4;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
    problem_data_0.omega_patches = 1:4;
    problem_data_0.gamma_sides = cell(9,1); 
    problem_data_0.gamma_sides(1:4) = {3, 3, 3, 3}; % relative to each patch
end