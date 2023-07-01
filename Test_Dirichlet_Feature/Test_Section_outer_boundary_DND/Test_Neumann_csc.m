clc
clearvars
close all
epsilon_ref = 1e-2 ;  
eps_values = epsilon_ref./ 2.^(0:6);
theta = @( x , y ) atan2( y , x );
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) cos( theta( x - 0.5 , y - 1 ) );%40*cos(pi*x).*cos(pi*y)+10*cos(5*pi*x).*cos(7*pi*y);
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 0 * x .* y ;
[problem_data, problem_data_0] = determineBC(problem_data);

ext = 4;
plotIt = 'true' ; 



method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [ 16 16 ]; % [32 32]; 
method_data.nquad = [5 5];


nn = numel(eps_values);
errh1s = zeros(1,nn);
est = zeros(1,nn);

for ii = 1:nn
    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
    %leave the two Dirichlet data on the trasmission regions here since they depend on epsilon
    
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
    [omega_0, msh_0, space_0, u_0] =  mp_solve_laplace_defeaturing_(problem_data_0, method_data);
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
%     normu(ii) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
%     errh1s_rel(ii) = errh1s(ii)/normu(ii);
%  
    %set the right Dirichlet data on \gamma: h_2 
    [est(ii), meas_gamma, err_l2_gamma , avg_d_gamma] = estimator_negative_dnd(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data.g,...
                        problem_data_0.omega_patches, problem_data.gamma_sides, ...
                        problem_data.omega0_patches, msh, space, u );
    
end


close all 

%%
if plotIt
    fig = figure;
    loglog(eps_values, errh1s, '+-r' , 'LineWidth', 4);
    hold on 
    loglog(  eps_values, est, '+-b' , 'LineWidth', 4)
    loglog( eps_values, eps_values, 'k:' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$', 'Estimator', '$\epsilon$'  , 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',35)
    title(['$g_{d}=0\:' ...
        '\gamma_0,\:\:\:\:g_{n}\:\gamma\cup\Gamma$'] , 'Interpreter','latex' , 'FontSize', 30 ...
    )
%     if saveIt
%         saveas(fig, filename, 'epsc');
%     end
end











function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    s = epsilon/(2*sqrt(2));
    srf_0(1) = nrbdegelev(nrbsquare([0.5-s, 1-s], 2*s, s), [1,1]);
    srf_0(1) = nrbkntins(srf_0(1), {0.5, 0.5});
    
    crv1 = nrbline([0.5-s, 1], [0.5-s, 1-s]);
    crv2 = nrbcirc(epsilon, [0.5, 1], pi, 5*pi/4);
    srf_0(2) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, 0.5});
    
    crv1 = nrbline([0.5-s, 1-s], [0.5+s, 1-s]);
    crv2 = nrbcirc(epsilon, [0.5, 1], 5*pi/4, 7*pi/4);
    srf_0(3) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(3) = nrbkntins(srf_0(3), {0.5, 0.5});
    
    crv1 = nrbline([0.5+s, 1-s], [0.5+s, 1]);
    crv2 = nrbcirc(epsilon, [0.5, 1], -pi/4, 0);
    srf_0(4) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(4) = nrbkntins(srf_0(4), {0.5, 0.5});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], pi, 5*pi/4);
    crv2 = nrbcirc( 1  , [0.5, 1], pi, 5*pi/4);
    srf_0(5) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(5) = nrbkntins(srf_0(5), {0.5, (ext-1)*epsilon/(1-epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], 5*pi/4, 7*pi/4);
    crv2 = nrbcirc( 1  , [0.5, 1], 5*pi/4, 7*pi/4);
    srf_0(6) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(6) = nrbkntins(srf_0(6), {0.5, (ext-1)*epsilon/(1-epsilon)});
    
    crv1 = nrbcirc(epsilon, [0.5, 1], -pi/4, 0);
    crv2 = nrbcirc( 1  , [0.5, 1], -pi/4, 0 );
    srf_0(7) = nrbdegelev(nrbruled(crv1, crv2), [0,1]);
    srf_0(7) = nrbkntins(srf_0(7), {0.5, (ext-1)*epsilon/(1-epsilon)});
%     for ii = 1 : 7
%         nrbplot( srf_0(ii ) , [ 20 20 ] )
%         hold on 
%         view( 2 )
%     end

    srf_F = srf_0(1:4);
    srf = srf_0(5:7);

%     nrbexport( srf , strcat( 'complete_problem.txt' ) )
%     nrbexport( srf_0 , strcat( 'defeatured_problem.txt' ) )
end


function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [2 4 7 3 5 8 ];
    problem_data.drchlt_sides_homogeneous = [];
    problem_data.drchlt_sides_inhomogeneous = [];
    problem_data.drchlt_sides_inhomogeneous_trasmission = [1 6];
    problem_data.drchlt_sides = [ 1 6 ] ;
    problem_data.gamma_sides = [2 4 7];
    problem_data.omega0_patches = 1:4;
%     problem_data.gamma_info = {[2 3 4], {4, 4, 4}, [1 3; 2 3; 3 3]};

    % Simplified problem
    problem_data_0.nmnn_sides = [ 5 6 8 ];
    problem_data_0.drchlt_sides_homogeneous =  [  ] ;
    problem_data_0.drchlt_sides_inhomogeneous = [];
    problem_data_0.drchlt_sides_inhomogeneous_trasmission = [1 2 3 4 7 ];
    problem_data_0.drchlt_sides = [   1 2 3 4 7 ] ;
    problem_data_0.omega_patches = 5:7;
    problem_data_0.gamma_sides = cell(7,1); 
    problem_data_0.gamma_sides(5:7) = {3, 3, 3}; % relative to each patch
end