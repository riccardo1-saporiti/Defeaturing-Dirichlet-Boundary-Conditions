clear
close all
clc
epsilon_ref = 1e-2;  
eps_values =  epsilon_ref ./ 2.^(-5:3); % 0.4./2.^(0:3); % 1e-2 ./ 2.^(-5:6); % 0.4./2.^(0:3); %0.1;
filename = 'results/test19_drchlt0';
saveIt = false;
plotIt = true;
problem_data.c_diff = @(x, y) ones(size(x));
problem_data.g = @(x, y, ind) zeros(size(x));
r = @( x , y ) sqrt( x .^ 2 + y .^ 2 );
theta = @( x , y ) atan2( y , x );
problem_data.f = @(x, y) 2 * r( x , y ) .* cos( theta( x , y ) ); 
[problem_data, problem_data_0] = determineBC(problem_data);
problem_data.h = @(x, y, ind) 0 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y   ;
problem_data_0.h = @(x, y, ind) 0 * x .* y ;
problem_data_0.h_auxiliar = @(x, y, ind) 0 * ( isempty( find( ind == [problem_data_0.gamma_sides{:}] ) ) == 0  ) + 0 * x .* y  ;
perturbation_type = 'singular' ; 
ext = 4;

method_data.degree = [3 3];
method_data.regularity = [2 2];
method_data.nsub = [ 3 3 ];
method_data.nquad = [5 5];

%%
neps = numel(eps_values);
errh1s = zeros(1,neps);
est = zeros(1,neps);
est_h1_tg = zeros( 1 , neps ) ; 
ref_vals = [ 3 , 16 , 64 ] ; 

for j = 1 : length( ref_vals )
    method_data.nsub = [ ref_vals( j ) , ref_vals( j ) ];
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
    
        %   set the right Dirichlet data on \gamma: h_2 
         [est( ii ), est_h1_tg( ii )  ] = estimator_negative_dirichlet(msh_0, space_0, u_0, problem_data_0.gamma_sides, problem_data_0.h_auxiliar,...
                            problem_data_0.omega_patches, space , epsilon  , perturbation_type  );
                                      
    end
    str_lab = strcat( 'c_test_star_shaped_' , num2str( ref_vals( j ) ) , '.mat' ) ; 
    save( str_lab , 'errh1s' , 'est' , 'est_h1_tg' ) 

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
    title(['$g_{d}=0\:\gamma\:\:\:\:f=2rcos(\theta)\:\Omega\:\:\:\:Shape:star\:\:\:\:64\:subdivisions$'] , 'Interpreter','latex' , 'FontSize', 30 ...
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

function [problem_data, problem_data_0] = determineBC(problem_data)
    problem_data_0 = problem_data;
    nn = 10;
    
    % Exact problem
    problem_data.nmnn_sides = []; 
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