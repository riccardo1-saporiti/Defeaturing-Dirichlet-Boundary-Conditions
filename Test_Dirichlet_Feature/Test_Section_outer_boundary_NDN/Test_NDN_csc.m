clc
clearvars
close all
epsilon_ref = 1e-2 ;  
eps_values = epsilon_ref./ 2.^(0:6);
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind)  0 .* x .* y ;
r = @( x , y ) sqrt( x .^ 2 + y .^ 2 );
theta = @( x , y ) atan2( y , x );
problem_data.f = @(x, y)  + 0 * x.* y ; %if the forcing term is odd in the x direction then we have automatic satisfaction of dirichlet b.c. at the origin
theta = @( x , y ) atan2( y , x );
[problem_data, problem_data_0] = determineBC(problem_data);
%impose nonhomogeneous Dirichlet data over the semi-circle \gamma
problem_data.h = @(x, y, ind) exp( sin( theta( x - 0.5 ,  y - 1 ) ) ) * ( isempty( find( ind == problem_data.drchlt_sides_inhomogeneous ) ) == 0  ) + 0 * x .* y  ;
problem_data_0.h = @(x, y, ind) 0 * x .* y ;
problem_data_0.h_auxiliar = @( x,  y , ind ) exp( sin( theta( x - 0.5 ,  y - 1 ) ) ) * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;

perturbation_type = 'singular' ;
ext = 4;
plotIt = 'true' ; 
method_data.degree = [ 3 3 ];
method_data.regularity = [ 2 2 ];
method_data.nsub = [ 8 8 ]; 
method_data.nquad = [5 5];


nn = numel(eps_values);
errh1s = zeros(1,nn);
est_h1_tg = zeros(1,nn);
est = zeros(1,nn);


for ii = 1:nn
        epsilon = eps_values(ii);
        fprintf('----- epsilon = %f -----\n', epsilon);

        % 1) BUILD GEOMETRY
        [srf_0, srf, srf_F] = buildGeometry(epsilon, ext);
        
        problem_data.geo_name = srf;
        problem_data_0.geo_name = srf_0;
        
        % 2) SOLVE THE EXACT PROBLEM
        [omega, msh, space, u] = mp_solve_laplace_defeaturing_ (problem_data, method_data);
        
        % strm = struct( msh );
        % for iiii = 1:3
        %     vect_plot = struct( strm.msh_patch{iiii} ).breaks ;
        %     sp_plot_solution (zeros( space.sp_patch{ iiii }.ndof , 1 ), space.sp_patch{ iiii }, omega( iiii ), vect_plot)
        %     view( 2 )
        %     hold on
        % end

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
    %     normu(ii) = errh1s_negative(msh, space, u, msh_0, space_0, zeros(size(u_0)), problem_data_0.omega_patches);
    %     errh1s_rel(ii) = errh1s(ii)/normu(ii);
    %  

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
    loglog( eps_values , 1 ./ sqrt( abs( log(eps_values) ) ) , '+-' , 'LineWidth' , 5 )
    grid on
    legend('$|u-u_0|_{1,\Omega}$', 'Estimator', '$Estimator(\nabla_{t}u_0)$' ,'$\frac{1}{\sqrt{ln(|\epsilon|)}}$' , 'Interpreter' , 'Latex' , 'Location', 'northwest' , 'Fontsize' , 30 )
    xlabel( '$\epsilon$' , 'Interpreter','latex' , 'FontSize', 30 )
    set(gca,'FontSize',30)
    title(['$NDN\rightarrow NNN,\:Circular\:section\:\:\:\:g_{d}=cos( x )\:\gamma\:\:\:\:g_{n}=0\:\:\:\:g_d=0\:\Gamma_d\setminus\gamma\:\:\:f=exp(-|x|^2)\:\:\:x^{4}\:refinement\:\:\:$'] , 'Interpreter','latex' , 'FontSize', 30 ...
    )
%     if saveIt
%         saveas(fig, filename, 'epsc');
%     end
end


ord_err = log2( errh1s( 1 : end - 1 ) ./ errh1s( 2 : end ) )
ord_est = log2( est( 1 : end - 1 ) ./ est( 2 : end ) )












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
    problem_data.nmnn_sides = [1 6];
    problem_data.drchlt_sides_homogeneous = [3 5 8];
    problem_data.drchlt_sides_inhomogeneous = [2 4 7];
    problem_data.drchlt_sides = [ 2 4 7 3 5 8 ] ;
    problem_data.gamma_sides = [2 4 7];
    problem_data.omega0_patches = 1:3;
%     problem_data.gamma_info = {[2 3 4], {4, 4, 4}, [1 3; 2 3; 3 3]};

    % Simplified problem
    problem_data_0.nmnn_sides = [1 2 3 4 7];
    problem_data_0.drchlt_sides_homogeneous =  [ 5 6 8 ] ;
    problem_data_0.drchlt_sides_inhomogeneous = [];
    problem_data_0.drchlt_sides= [ 5 6 8 ];
    problem_data_0.omega_patches = 5:7;
    problem_data_0.gamma_sides = cell(7,1); 
    problem_data_0.gamma_sides(5:7) = {3, 3, 3}; % relative to each patch
end
