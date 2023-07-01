clc
clearvars
close all
epsilon_ref = 1e-2;
eps_values = epsilon_ref./ 2.^(0:6);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) 40*cos(pi*x).*cos(pi*y)+10*cos(5*pi*x).*cos(7*pi*y);
problem_data.f = @(x, y) zeros(size(x)); 
problem_data.h = @(x, y, ind) 0 * x .* y ;
theta = @( x , y ) atan2( y , x );
[problem_data, problem_data_0] = determineBC(problem_data);
ext = 4;
plotIt = 'true' ; 

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [  8 8 ]; % [32 32]; 
method_data.nquad = [5 5];


nn = numel(eps_values);
errh1s = zeros(1,nn);
errh1s_interface = zeros(1,nn);
est = zeros(1,nn);
normu = zeros(1,nn);
errh1s_rel = zeros(1,nn);
meas_gamma = zeros(1,nn);


for ii =  1 : nn

    epsilon = eps_values(ii);
    fprintf('----- epsilon = %f -----\n', epsilon);
   
    % 1) BUILD GEOMETRY
    [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
    problem_data.geo_name = srf;
    problem_data_0.geo_name = srf_0;
    
    % 2) SOLVE THE EXACT PROBLEM
    [omega, msh, space, u] = mp_solve_laplace_defeaturing_ (problem_data, method_data);
    %check the mesh 
    % strm = struct( msh ); 
    % for ii = 1:length( strm.msh_patch )
    %     vect_plot = struct( strm.msh_patch{ii} ).breaks ;
    %     sp_plot_solution (u * 0, space, omega, vect_plot)
    %     hold on
    % end
    vtk_pts = {linspace(0, 1, 50), linspace(0, 1, 50)};
    figure
    subplot( 1 , 2 , 1 )
    sp_plot_solution (u, space, omega, vtk_pts)
    % 3) SOLVE THE DEFEATURED PROBLEM
    [omega_0, msh_0, space_0, u_0] = mp_solve_laplace_defeaturing_ (problem_data_0, method_data);

    % strm = struct( msh_0 ); 
    % for ii = 1:length( strm.msh_patch )
    %     vect_plot = struct( strm.msh_patch{ii} ).breaks ;
    %     sp_plot_solution (u_0 * 0, space_0, omega_0, vect_plot)
    %     hold on
    % end

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
    loglog( eps_values , errh1s , '+-r' , 'LineWidth', 3);
    hold on 
    loglog(  eps_values, est , 'b+-' , 'LineWidth', 3)
    loglog( eps_values, eps_values , 'g:' , 'LineWidth' , 5 )
    loglog( eps_values, eps_values .^ 2, 'k:' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$' , '$Estimator(u_0)$' , '$\epsilon$' , '$\epsilon^2$', 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$g_{d}=0\:' ...
        '\partial{\Omega}\setminus\gamma,\:\:\:\:\:g_{n}\:\gamma\:\:Squared\:section\:square\:(central)$'] , 'Interpreter','latex' , 'FontSize', 50 ...
    )
%     if saveIt
%         saveas(fig, filename, 'epsc');
%     end
end





%%
function [srf_0, srf, srf_F] = buildGeometry(epsilon, ext)
    srf_0(1) = nrb4surf ([0 0], [0.5-epsilon/2 0], [0 1-epsilon], [0.5-epsilon/2 1-epsilon]);
    srf_0(1) = nrbkntins(srf_0(1), {(0.5-ext*epsilon/2)/(0.5-epsilon/2), (1-ext*epsilon/2)/(1-epsilon/2)});

    srf_0(2) = nrb4surf ([0.5-epsilon/2 0], [0.5+epsilon/2 0], [0.5-epsilon/2 1-epsilon], [0.5+epsilon/2 1-epsilon]);
    srf_0(2) = nrbkntins(srf_0(2), {0.5, (1-ext*epsilon/2)/(1-epsilon/2)});

    srf_0(3) = nrb4surf ([0.5+epsilon/2 0], [1 0], [0.5+epsilon/2 1-epsilon], [1 1-epsilon]);
    srf_0(3) = nrbkntins(srf_0(3), {1-(1-ext*epsilon)/(1-epsilon), (1-ext*epsilon/2)/(1-epsilon/2)});

    srf_0(4) = nrb4surf ([0 1-epsilon], [0.5-epsilon/2 1-epsilon], [0 1], [0.5-epsilon/2 1]);
    srf_0(4) = nrbkntins(srf_0(4), {(0.5-ext*epsilon/2)/(0.5-epsilon/2), 0.5});

    srf_0(5) = nrb4surf ([0.5-epsilon/2 1-epsilon], [0.5+epsilon/2 1-epsilon], [0.5-epsilon/2 1], [0.5+epsilon/2 1]);
    srf_0(5) = nrbkntins(srf_0(5), {0.5, 0.5});

    srf_0(6) = nrb4surf ([0.5+epsilon/2 1-epsilon], [1 1-epsilon], [0.5+epsilon/2 1], [1 1]);
    srf_0(6) = nrbkntins(srf_0(6), {1-(1-ext*epsilon)/(1-epsilon), 0.5});

    srf = srf_0([1:4, 6]);
    srf_F = srf_0(5);
    
    for ii = 1 : 5
        nrbplot( srf( ii ) , [20 20] );
        hold on
        view( 2 )
    end


end

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    
    % Exact problem
    problem_data.nmnn_sides = [ 4 8 10 ];
    problem_data.drchlt_sides = [ 7 1 2 3 6 5 11 9 12 ];
    problem_data.gamma_sides = [4 8 10];
    problem_data.omega0_patches = 1:5;

    % Simplified problem
    problem_data_0.nmnn_sides = [ ];
    problem_data_0.drchlt_sides = [ 6 1 2 3 5 4 9  8 7 10]  ; 
    problem_data_0.omega_patches = [1:4 6];
    problem_data_0.gamma_sides = cell(6,1); 
    problem_data_0.gamma_sides([2 4 6]) = {4, 2, 1}; % relative to each patch
    problem_data_0.patches_gamma_0 = 5 ;
    problem_data_0.gamma_0_sides = cell( 6 , 1 );
    problem_data_0.gamma_0_sides( 5 ) = { 4 };
end
