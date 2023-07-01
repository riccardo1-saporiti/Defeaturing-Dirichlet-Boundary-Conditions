clc
clearvars
close all
epsilon_ref = 1e-2;  
eps_values = epsilon_ref./ 2.^(0:6);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y)  +  0 * x .* y ; 
theta = @( x , y ) atan2( y , x );
r = @( x , y ) sqrt( x .^ 2 + y .^ 2 );
[problem_data, problem_data_0] = determineBC(problem_data);
problem_data.h = @(x, y, ind) 1 * ( isempty( find( ind == problem_data.drchlt_sides_inhomogeneous ) ) == 0  ) + 0 * x .* y  ;
problem_data_0.h = @(x, y, ind) 0 * x .* y  ;
problem_data_0.h_auxiliar = @( x , y , ind ) 1 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;

ext = 4;
plotIt = 'true' ; 


method_data.degree = [3 3];
method_data.regularity = [2 2 ];
method_data.nsub = [ 8 8 ]; % [32 32]; 
method_data.nquad = [5 5];
problem_data_0.p = method_data.degree(  1 ) ;  

nn = numel(eps_values);
errh1s = zeros(1,nn);
errh1s_interface = zeros(1,nn);
est = zeros(1,nn);
est_h1_tg = est ; 

perturbation_type = 'singular' ; 
problem_data_0.singularity_type = 'const' ;
problem_data_0.dirac_center = [ 0 ; 0 ; 0 ];
problem_data_0.Dirac_Test = 'true';
Find_alpha_opt_Dirac = 1;

for ii = nn:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    problem_data_0.epsilon = epsilon ; 
    problem_data_0.alpha_dirac =  - 2 * pi / log( epsilon );
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);


    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;

    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_defeaturing_ (problem_data, method_data);
    
    % % 
    % % % 
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
    figure
    subplot( 1 , 2 , 1 )
    sp_plot_solution (u, space, omega, vtk_pts)
    
    if Find_alpha_opt_Dirac
        problem_data.u_0_exact = @( x , y , alpha ) -alpha * log( sqrt( x .^2 + y.^ 2 ) ) / ( 2*pi) ;
        [~, boundaries, ~, ~, ~] = mp_geo_load (problem_data.geo_name);
        gamma_side= problem_data.gamma_sides ; 
        l2_norm_bdry = @( alpha ) get_l2_norm_boundary(boundaries,space,msh,gamma_side,problem_data,alpha);
        alpha_star = fsolve( l2_norm_bdry , problem_data_0.alpha_dirac) ;
        problem_data_0.alpha_dirac = alpha_star; 
    end
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_defeaturing_ (problem_data_0, method_data);
    %modify this function to consider the new mesh for the point 0,0,0
    % strm = struct( msh_0 );
    % for ii = 1:7
    %     vect_plot = struct( strm.msh_patch{ii} ).breaks ;
    %     sp_plot_solution (u_0 * 0, space_0, omega_0, vect_plot)
    % 
    % end
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
    subplot( 1 , 2 , 2 )
    sp_plot_solution ( u_0 , space_0 , omega_0 , vtk_pts)
    
    % check how good is the approximation of the Delta computing the l2 and
    % semi h1 error of the exact defeatured solution and its approximation
    % syms rs cos_th x y
    % rs = sqrt( x ^ 2 + y ^ 2 ) ;
    % cos_th = x / rs;
    % f_ex =  - problem_data_0.alpha_dirac * log( rs ) / ( 2 * pi );
    % df_x = simplify( diff( f_ex , x ) );
    % df_y = simplify( diff( f_ex , y ) );
    % problem_data.uex_def = matlabFunction( f_ex);  %exsolution
    % problem_data.graduex_def = @(x,y ) ( cat( 1 , reshape( eval( df_x ) , [ 1 , size( x ) ] ) , reshape( eval( df_y ) , [ 1 , size( x ) ] ) ) );
    % [ ~ , errl2_def_var , errh1s_def_var ] = sp_h1_error (space_0, msh_0, u_0, problem_data.uex_def, problem_data.graduex_def );
    % [ ~ , errl2_uex_def_var , errh1s_uex_def_var ] = sp_h1_error (space, msh, u, problem_data.uex_def, problem_data.graduex_def );
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
    legend('$|u-u_0|_{1,\Omega}$', 'Estimator', '$Estimator(\nabla_{t}u_0)$'  , 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$g_{d}=1\:\gamma\:\:\:f=0\:\:' ...
        '\Omega$'] , 'Interpreter','latex' , 'FontSize', 30 ...
    )
%     if saveIt
%         saveas(fig, filename, 'epsc');
%     end
end

% save('dirac_test_32_refs_1e-5.mat' , 'errh1s' , 'est_h1_tg')

ord_eps_est_h1 = log2( est_h1_tg( 1 : end - 1 ) ./ est_h1_tg( 2 : end ) )
ord_eps_est = log2( est( 1 : end - 1 ) ./ est( 2 : end ) )
ord_eps_err = log2( errh1s( 1 : end - 1 ) ./ errh1s( 2 : end ) )
1 ./ sqrt( abs( log( eps_values) ) )

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
    % for ii = 1 : 8
    %     nrbplot( srf( ii ) , [ 20 20 ] );
    %     hold on
    % end
end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [];
    problem_data.drchlt_sides_homogeneous = 1:4;
    problem_data.drchlt_sides_inhomogeneous = 5:8;
    problem_data.drchlt_sides =  1:8 ;
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