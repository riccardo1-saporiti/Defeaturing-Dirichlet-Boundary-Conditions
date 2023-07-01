function [est , est_h1_tg  ] = estimator_negative_dirichlet(msh_0, space_0, u_0,...
            gamma_sides_local_forOmega0, gfun, inter_patches_forOmega0, space, epsilon , perturbation_type, lambda, mu )
    
    %Compute stiffness matrix with the tangential gradient and mass matrix
    A_gradvt_gradut = spalloc (space_0.ndof, space_0.ndof, 3*space_0.ndof);
    A_uv = spalloc (space_0.ndof, space_0.ndof, 3*space_0.ndof);
    gD = zeros(space_0.ndof, 1);
    cs = cumsum([0, msh_0.nel_per_patch]);
    internal_element_list_0 = [];
    for iptc = inter_patches_forOmega0 
        space_ptc = space_0.sp_patch{iptc};
        msh_ptc = msh_0.msh_patch{iptc};
        dofs_ptc = space_0.gnum{iptc};
        internal_element_list_0 = [internal_element_list_0 cs(iptc)+1:cs(iptc+1)];
        msh_inner = msh_evaluate_element_list (msh_0, internal_element_list_0);
        space_inner = sp_evaluate_element_list (space_0, msh_inner, 'value', 1, 'gradient', 1);
        for iside = gamma_sides_local_forOmega0{iptc}
            msh_side = msh_eval_boundary_side (msh_ptc, iside);
            msh_side_from_interior = msh_boundary_side_from_interior (msh_ptc, iside);

            sp_bnd = space_ptc.constructor (msh_side_from_interior);
            if space.ncomp == 1
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true);
                
            else
                sp_bnd = sp_precompute (sp_bnd, msh_side_from_interior, 'value', true, 'gradient', true, 'divergence', true);
            end

            coeff_at_qnodes = ones (msh_side.nqn, msh_side.nel);

            B_uv = op_u_v (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            A_uv(dofs_ptc,dofs_ptc) = A_uv(dofs_ptc,dofs_ptc) + B_uv;
            if space.ncomp == 1    
                B_gradvt_gradut = op_gradu_t_gradv_t (sp_bnd, sp_bnd, msh_side, coeff_at_qnodes);
            end
            A_gradvt_gradut(dofs_ptc,dofs_ptc) = A_gradvt_gradut(dofs_ptc,dofs_ptc) + B_gradvt_gradut;
            [gD_bddofs, bddofs] = sp_drchlt_l2_proj ( space_ptc , msh_ptc , gfun , iside );
            gD(dofs_ptc(bddofs)) = gD_bddofs;
        end

    end

    
   
    
 
    meas_gamma = ones(size(A_uv,2),1).' * A_uv * ones(size(A_uv,2),1) / space.ncomp;
 
    est_interface_avg2 = u_0.' * A_uv * u_0 ...
                          + gD.' * A_uv * gD ...
                          - 2*u_0.' * A_uv * gD;
    est_interface_avg2_semi_h1_tg = u_0' * A_gradvt_gradut * u_0 + gD.' * A_gradvt_gradut * gD - ...
        2 * u_0.' * A_gradvt_gradut * gD ; 
    
    if space.ncomp == 1
        integral_compat = gD.' * A_uv * ones(size(gD)) ...
            - u_0.' * A_uv * ones(size(gD));
    else
        dofs_per_comp = [];
        for icomp = 1:space_0.ncomp
            global_dofs = [];
            for iptc = 1:space_0.npatch
                local_dofs = space_0.sp_patch{iptc}.comp_dofs{icomp};
                global_dofs = [global_dofs space_0.gnum{iptc}(local_dofs)];
            end
            global_dofs = unique(global_dofs);
            dofs_per_comp = [dofs_per_comp; global_dofs];
        end

        integral_compat2 = 0;
        for icomp = 1:space_0.ncomp
            ei = zeros(size(gD));
            ei(dofs_per_comp(icomp,:)) = ones(size(dofs_per_comp(icomp,:)));
            integral_compat2 = integral_compat2 + ( gD.' * A_uv * ei ...
            - u_0.' * A_uv * ei )^2;
        end
        integral_compat = sqrt(abs( integral_compat2) );
    end
    est_interface2 = est_interface_avg2 - integral_compat^2 / meas_gamma;   %can check that this holds
    
    % To deal with round-off accuracy issues
    if isequal( (integral_compat^2 / meas_gamma) - [0 1e-10] > est_interface_avg2, [1 0])
        est_interface2 = abs(est_interface2);
    end
    
    switch perturbation_type
        case{'singular'}
            c_epsilon = ( 1 / sqrt( abs( log( abs( epsilon ) ) ) ) ) ;
        case{'regular'}
            c_epsilon = 1;
    end

    rho_epsilon = 1 / sqrt( epsilon ) ; 

    if msh_0.rdim == 2
        
        est = c_epsilon * abs( integral_compat / meas_gamma ) + rho_epsilon * sqrt( est_interface2 );
        est_h1_tg = c_epsilon * abs( integral_compat / meas_gamma ) + sqrt( sqrt( abs( est_interface2 * est_interface_avg2_semi_h1_tg ) ) );
        
    elseif msh_0.rdim == 3  %three dimensional case not yet studied
        
        est = c_epsilon * abs( integral_compat / meas_gamma ) + rho_epsilon * sqrt( est_interface2 );
        est_h1_tg = c_epsilon * abs( integral_compat / meas_gamma ) + sqrt( sqrt( abs( est_interface2 * est_interface_avg2_semi_h1_tg ) ) );
        
    else
        error('rdim should be 2 or 3')
    end
end


