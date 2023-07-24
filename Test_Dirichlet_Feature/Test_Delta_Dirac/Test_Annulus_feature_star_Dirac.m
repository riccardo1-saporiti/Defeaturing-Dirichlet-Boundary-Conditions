clc
clearvars
close all
epsilon_ref = 1e-2;  
eps_values =  epsilon_ref ./ 2.^(-5:3); % 0.4./2.^(0:3); % 1e-2 ./ 2.^(-5:6); % 0.4./2.^(0:3); %0.1;
filename = 'results/test19_drchlt0';
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.f = @(x, y)  +  0 * x .* y ; 
theta = @( x , y ) atan2( y , x );
r = @( x , y ) sqrt( x .^ 2 + y .^ 2 );
[problem_data, problem_data_0 , problem_data_F] = determineBC(problem_data);
problem_data.h = @(x, y, ind) 1 * ( isempty( find( ind == problem_data.drchlt_sides_inhomogeneous ) ) == 0  ) + 0 * x .* y  ;
problem_data_0.h = @(x, y, ind) 0 * x .* y  ;
problem_data_0.h_auxiliar = @( x , y , ind ) 1 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;

ext = 4;
plotIt = 'true' ; 


method_data.degree = [3 3];
method_data.regularity = [2 2 ];
method_data.nsub = [ 3 3 ]; % [32 32]; 
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
J_alpha = 'H^{1/2}';

for ii = 1:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    
    problem_data_0.epsilon = epsilon ; 
    problem_data_0.alpha_dirac =  - 2 * pi / log( epsilon );
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);


    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    problem_data_F.geo_name = srf_F;
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
%         alpha_star = fsolve( l2_norm_bdry , problem_data_0.alpha_dirac) ;
        switch J_alpha
            case{'L^{2}'}

                J_cost = @( alpha ) l2_norm_bdry( alpha );

            case{'H^{1/2}'}


                semi_h12_boundary_norm = @( alpha ) get_semi_h12_norm_boundary(problem_data_F , method_data , alpha );
%                 meas_gamma = 8 * epsilon;
                beta = 1;
                xi = 0.05;
                avg_term = @(alpha ) get_average_boundary(boundaries,space,msh,gamma_side,problem_data,meas_gamma,alpha);
                J_cost = @( alpha ) l2_norm_bdry( alpha ) + beta * semi_h12_boundary_norm( alpha ) ; 
        end
        alpha_star = fsolve( J_cost , problem_data_0.alpha_dirac) ;
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
    [est( ii ), est_h1_tg( ii ) , avg_d_gamma , h12est ] = estimator_negative_dirichlet(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data_0.h_auxiliar ,...
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
    perim = 4*epsilon;
%     area = epsilon^2;
    nn = 10;
    
%     rin = sqrt( area / (nn/2*sin(2*pi/nn)*( 1+(2-sin((nn-2)*pi/(2*nn)))/(sin(8*pi/20)) )) );
    rin = perim / (2*nn*sqrt(5-4*cos(2*pi/(2*nn))));
    rout = 2*rin;
    
    for ii = 0:nn-1
        thetain = ii * 2*pi/nn;
        if ii > 0
            line = nrbline([rout*cos(thetaout) rout*sin(thetaout)], [rin*cos(thetain) rin*sin(thetain)]);
            circ = nrbcirc(1, [0,0], thetaout, thetain);
            srf(2*ii) = nrbruled(line, circ);

            srf_F(ii) = nrbruled(nrbline([rin*cos((ii-1) * 2*pi/nn), rin*sin((ii-1) * 2*pi/nn)], [0,0]), line);
        end
        thetaout = thetain + pi/nn;
        line = nrbline([rin*cos(thetain) rin*sin(thetain)], [rout*cos(thetaout) rout*sin(thetaout)]);
        circ = nrbcirc(1, [0,0], thetain, thetaout);
        srf(2*ii+1) = nrbruled(line, circ);
    end
    thetain = 0;
    line = nrbline([rout*cos(thetaout) rout*sin(thetaout)], [rin*cos(thetain) rin*sin(thetain)]);
    circ = nrbcirc(1, [0,0], thetaout, thetain);
    srf(2*nn) = nrbruled(line, circ);
    
    srf_F(nn) = nrbruled(nrbline([rin*cos((nn-1) * 2*pi/nn), rin*sin((nn-1) * 2*pi/nn)], [0,0]), line);
    for ii = 1:nn
        srf_F(ii) = nrbdegelev(srf_F(ii), [1, 1]);
%         ktF = max(0.5, 1-ext);
        srf_F(ii) = nrbkntins(srf_F(ii), {0.5, 0.5});
        srf(2*ii-1) = nrbdegelev(srf(2*ii-1), [0, 1]);
        srf(2*ii) = nrbdegelev(srf(2*ii), [0, 1]);
%         kt = min(0.5, 3/2*rin*ext/(1-3/2*rin));
        kt = linspace(0,1,10);
        srf(2*ii-1) = nrbkntins(srf(2*ii-1), {0.5, kt});
        srf(2*ii) = nrbkntins(srf(2*ii), {0.5, kt});
    end

    srf_0 = srf;
    srf_0(2*nn+(1:nn)) = srf_F;
    % for ii = [ 1:20 ]
    %     nrbplot( srf(ii ) , [ 20 20 20 ] )
    %         hold on 
    % 
    % end
%     ns = 2;
%     nsub = cell(3*nn, 1);
%     [nsub{1:2*nn}] = deal([ns, ns*10]); 
%     [nsub{2*nn+1:end}] = deal([ns, ns]);
end

function [problem_data, problem_data_0 , problem_data_F] = determineBC(problem_data)
    problem_data_0 = problem_data;
    problem_data_F = problem_data;
    problem_data_F.nmnn_sides = [];
    problem_data_F.drchlt_sides = 1:20;

    nn = 10;
    
    % Exact problem
    problem_data.nmnn_sides = []; 
    problem_data.drchlt_sides_homogeneous = 2:2:4*nn;
    problem_data.drchlt_sides_inhomogeneous = 1:2:4*nn;
    problem_data.drchlt_sides = [ 1:2:4*nn  2:2:4*nn ];
    problem_data.gamma_sides = 1:2:4*nn;
    problem_data.omega0_patches = 1:2*nn;

    % Simplified problem
    problem_data_0.nmnn_sides = [];
    problem_data_0.drchlt_sides = 1:2*nn;
    problem_data_0.omega_patches = 1:2*nn;
    problem_data_0.gamma_sides = cell(3*nn,1); 
    [problem_data_0.gamma_sides{1:2*nn}] = deal(3); % relative to each patch
end

