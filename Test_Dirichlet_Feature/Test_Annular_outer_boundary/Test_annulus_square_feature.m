clear
close all
clc
epsilon_ref = 1e-2 ; 
eps_values = epsilon_ref ./ 2.^(0:6); 
filename = 'results/test18_drchlt0';
saveIt = false;
plotIt = true;
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
r = @( x , y ) sqrt( x .^ 2 + y .^ 2 );
theta = @( x , y ) atan2( y , x );
problem_data.f = @(x, y) 2 * r( x , y ) .* cos( theta( x , y ) ); 
problem_data.h = @(x, y, ind) 0 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y   ;
[problem_data, problem_data_0] = determineBC(problem_data); 
problem_data_0.h_auxiliar = @(x, y, ind) 0 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;
perturbation_type = 'singular' ;
ext = 4;
method_data.degree = [ 3 3 ];
method_data.regularity = [ 8 8];
method_data.nsub = [ 8 8 ];
method_data.nquad = [ 5 5 ];

%%
neps = numel(eps_values);
errh1s = zeros(1,neps);
errh1s_interface = zeros(1,neps);
est = zeros(1,neps);
est_h1_tg = est ; 

for ii = 1:neps
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
    subplot( 1 , 2 , 2 )
    sp_plot_solution (u_0 , space_0 , omega_0 , vtk_pts)

    % 4) COMPUTE ERROR AND ESTIMATOR
     errh1s(ii) = errh1s_negative(msh, space, u, msh_0, space_0, u_0, problem_data_0.omega_patches);

    %set the right Dirichlet data on \gamma: h_2 
     [est( ii ), est_h1_tg( ii )  ] = estimator_negative_dirichlet(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data_0.h_auxiliar,...
                        problem_data_0.omega_patches, space , epsilon , perturbation_type  );
                       
end

%%
close all
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r' , 'LineWidth', 4);
    hold on 
    loglog(  eps_values, est, '+-b' ,  'LineWidth', 4)
    loglog( eps_values , est_h1_tg , '+-' , 'LineWidth', 4  )
    loglog( eps_values, eps_values, 'k:' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$', 'Estimator', '$Estimator(\nabla_{t}u_0)$' , '$\epsilon$' , 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$g_{d}=0\:\gamma\:\:\:\:f=2rcos(\theta)\:\Omega\:\:\:\:Shape:square$'] , 'Interpreter','latex' , 'FontSize', 30 ...
    )
end


%%

function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    L = epsilon/2;
    tol = 0.1;
    srf_0(1) = nrbruled (nrbcirc(tol, [0 0], pi/4, 3*pi/4), nrbcirc(1, [0 0], pi/4, 3*pi/4));
    srf_0(2) = nrbtform(srf_0(1), vecrotz(pi/2));
    srf_0(3) = nrbtform(srf_0(2), vecrotz(pi/2));
    srf_0(4) = nrbtform(srf_0(3), vecrotz(pi/2));
    
    srf_0(5) = nrbruled (nrbline([L L], [-L L]), nrbcirc(tol, [0 0], pi/4, 3*pi/4));
    srf_0(6) = nrbtform(srf_0(5), vecrotz(pi/2));
    srf_0(7) = nrbtform(srf_0(6), vecrotz(pi/2));
    srf_0(8) = nrbtform(srf_0(7), vecrotz(pi/2));

    if ext*epsilon < 0.1
        for ii = 1:4
            srf_0(ii) = nrbkntins(srf_0(ii), {0.5, 0.5});
            srf_0(ii+4) = nrbkntins(srf_0(ii+4), {0.5, ext*epsilon/0.1});
        end
    else
        for ii = 1:4
            srf_0(ii) = nrbkntins(srf_0(ii), {0.5, (ext*epsilon-0.1)/0.9});
            srf_0(ii+4) = nrbkntins(srf_0(ii+4), {0.5, 0.5});
        end
    end
    srf_0(9) = nrbdegelev(nrb4surf([-L -L], [L, -L], [-L L], [L L]), [1 1]);
    srf_0(9) = nrbkntins(srf_0(9), {0.5, 0.5});
    
    srf = srf_0(1:8);
    srf_F = srf_0(9);
    % for ii = 5 : 8
    %     nrbplot( srf_0( ii ) , [ 20 20 ] );
    %     hold on
    %     view ( 2)
    % end
    % title( 'Patches 5:8 of $\Omega$' , 'Interpreter','latex','FontSize', 35)
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [];
    problem_data.drchlt_sides = 1:8;
    problem_data.gamma_sides = 5:8;
    problem_data.omega0_patches = 1:8;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:4;
    problem_data_0.omega_patches = 1:8;
    problem_data_0.gamma_sides = cell(9,1); 
    problem_data_0.gamma_sides(5:8) = {3, 3, 3, 3}; % relative to each patch
end

function res = myf(x, y)
    tol = 0.1;
    r = sqrt(x.^2+y.^2);
    res = ones(size(x));
    res(r<=tol) = zeros(size(res(r<=tol)));
end