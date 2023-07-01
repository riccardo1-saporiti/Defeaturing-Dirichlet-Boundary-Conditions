% SP_TO_VTK: Export to VTK format for plotting.
%
%  sp_to_vtk (u, space, msh, geometry, npts_per_el, filename)
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object representing the space of discrete functions (see sp_trimmed)
%     msh  :       object representing the trimmed mesh (see msh_trimmed)
%     geometry:    geometry structure (see geo_load)
%     npts_per_el: number of evaluation points on each element, along each parametric direction
%     filename:    name of the output file. 
%
% OUTPUT:
%
%    none
%
% Warning: the use of the function is different from the one in the class
%           sp_scalar.
% 
% Copyright (C) 2020, 2021 Luca Coradello
% Copyright (C) 2021 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function sp_to_vtk (u, sp_trimmed, msh_trimmed, geometry, npts_per_el, filename)

F = [];
solution = [];
active_el_index = 0; 
geo_tolerance = 1e-8;

trimmed_elems = msh_trimmed.reparam.trimmed_elems;    
non_trimmed_ids = msh_trimmed.reparam.non_trimmed_elem_ids;  

ndim = msh_trimmed.ndim;
global_to_active = sp_trimmed.global_to_active;

nel_dir = msh_trimmed.msh_cart.nel_dir;

for i_untrimmed_el = 1:numel(non_trimmed_ids)

        [I,J] = ind2sub(nel_dir, non_trimmed_ids(i_untrimmed_el)); 

        linear_ind = [I,J];

        knt = cell (ndim, 1);
        for idim=1:ndim
            knt{idim} = msh_trimmed.msh_cart.breaks{idim}(linear_ind(idim):linear_ind(idim)+1);
        end

        for idim = 1:ndim
            pts{idim} = linspace (knt{idim}(1)+geo_tolerance, knt{idim}(end)-geo_tolerance, npts_per_el(idim))';
        end

        msh_dummy = msh_cartesian (knt, pts, [], geometry, 'boundary', false);
        sp  = sp_trimmed.space_untrimmed.constructor (msh_dummy);

        if(size(msh_dummy.nel_dir,2) == 1)
            msh_dummy.nel_dir = msh_dummy.nel_dir';
        end

        msh_struct = msh_precompute(msh_dummy);
        sp_struct = sp_precompute (sp, msh_dummy, 'value', true);

        u_el = u(global_to_active(sp_struct.connectivity));

        u_el = reshape(u_el, 1, []);

        sol_el = reshape(squeeze(sum(bsxfun(@times, sp_struct.shape_functions, u_el), 2 )), 1, []);

        F_el = msh_struct.geo_map; 

        F_el = reshape (F_el, [msh_trimmed.rdim, npts_per_el]);

        sol_el = reshape (sol_el, [1, npts_per_el]);

        F = cat(4, F, F_el);
        solution = cat(4, solution, sol_el);

        active_el_index = active_el_index + 1;

end


% trimmed_elems_of_level = [];
% if(~isempty(trimmed_elems_ids))
%     for iindex = 1:numel(IA)
%         trimmed_elems_of_level = [trimmed_elems_of_level trimmed_elems(IA(iindex))];
%     end
% end

for i_trimmed_el = 1:msh_trimmed.reparam.nb_trimmed_elems

        for jtile = 1:trimmed_elems(i_trimmed_el).nb_tiles

            srf = trimmed_elems(i_trimmed_el).tiles_geometry(jtile); % NURBS geometry of the patch

            zeta = {unique(srf.knots{1}), unique(srf.knots{2})};

            for idim = 1:ndim
                pts{idim} = linspace (zeta{idim}(1)+geo_tolerance, zeta{idim}(end)-geo_tolerance, npts_per_el(idim))';
            end

            for jj = 1:ndim
                brk{jj} = [zeta{jj}(1) zeta{jj}(end)];
            end

            mesh_tile = msh_cartesian (brk, pts, [], geo_load (srf), 'boundary', false); % msh_cartesian structure

%             mesh_tile.nel_dir = mesh_tile.nel_dir(:);

            msh_tile = msh_precompute(mesh_tile);

            quad_pnts_in_srf_sp = msh_tile.geo_map;

            F_tile = zeros (msh_trimmed.rdim, msh_tile.nqn);

            for iPoint = 1:size(quad_pnts_in_srf_sp, 2)
                aux_brk{1} = [0 1];
                aux_brk{2} = [0 1];

                aux_point{1} = quad_pnts_in_srf_sp(1,iPoint);
                aux_point{2} = quad_pnts_in_srf_sp(2,iPoint);

                mesh_dummy = msh_cartesian (aux_brk, aux_point, [], geometry, 'boundary', false); % msh_cartesian structure
                sp  = sp_trimmed.space_untrimmed.constructor (mesh_dummy);

                msh = msh_precompute(mesh_dummy);
                sp_struct = sp_precompute (sp, msh, 'value', true);

                F_tile(:,iPoint) = msh.geo_map;

                u_tile = u(global_to_active(sp_struct.connectivity));

                u_tile = reshape(u_tile, 1, []);

                sol_tile(:,iPoint) = reshape(squeeze(sum(bsxfun(@times, sp_struct.shape_functions, u_tile), 2 )), 1, []);

            end

            F_tile = reshape (F_tile, [msh_trimmed.rdim, npts_per_el]);
            sol_tile = reshape (sol_tile, [1, npts_per_el]);

            F = cat(4, F, F_tile);

            solution = cat(4, solution, sol_tile);

            active_el_index = active_el_index + 1;

            clear F_tile bending_tile membrane_tile

        end


end

data_title = 'solution plot';

if(msh_trimmed.rdim == 2)
    x = reshape(F(1,:,:,:),[],1);
    y = reshape(F(2,:,:,:),[],1);
    z = zeros(size(x));
elseif(msh_trimmed.rdim == 3)
    x = reshape(F(1,:,:,:),[],1);
    y = reshape(F(2,:,:,:),[],1);
    z = reshape(F(3,:,:,:),[],1);   
else
    error('not implemented')
end

data_struct(1).name = 'solution';
data_struct(1).type = 'scalar';

u_1 = reshape(solution(1,:,:,:),[],1);

data_struct(1).data = [u_1]';

nel_dir = npts_per_el - 1;

cells_per_el = prod(nel_dir);
connectivity_el = zeros(4, cells_per_el);
index = 1;

npts_el = prod(npts_per_el);

for j = 1:nel_dir(2)
    
    m = npts_per_el(1);
    
    for i = 1:nel_dir(1)
        
        connectivity_el(1,index) = i + m*(j-1);
        connectivity_el(2,index) = i+1 + m*(j-1);
        connectivity_el(3,index) = i + m*j;
        connectivity_el(4,index) = i+1 + m*j;
        
        index = index + 1;
    end
    
end

connectivity = zeros(4,cells_per_el*active_el_index);
index = 1; shifting = 0;
for i = 1:active_el_index
    connectivity(:,index:(index+cells_per_el-1)) = connectivity_el + shifting;
    index = index + cells_per_el; 
    shifting = shifting + npts_el;
end

data_points = [x y z]';

vtk_write_quad_grid_and_data(filename, data_title, data_points, connectivity, data_struct, true);

end

