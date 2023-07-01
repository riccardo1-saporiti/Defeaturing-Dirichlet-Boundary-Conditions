% SP_EVALUATE_TILES_PARAM: create a struct of a scalar space defined over a tile coming from the reparametrization of a trimmed element,
%                                       where the evaluation is performed in the parametric domain
% 
%     sp_tile = sp_evaluate_tiles_param (space, msh_tile, varargin)
%     
% INPUTS:
%     
%   space:     object representing the space (see sp_trimming)
%   msh_tiles: struct describing the mesh info associated to tiles, see msh_evaluate_tiles
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%
% OUTPUT:
%
%    sp_tile: struct representing the discrete function space, with the
%    fields described in sp_evaluate_element_list_param
%              (see the article for a detailed description)
%
% Copyright (C) 2020 Ondine Chanon, Luca Coradello, Riccardo Puppi
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

function  sp_tile = sp_evaluate_tiles_param (space, msh_tile, varargin)

value = true;
gradient = false;
hessian = false;

if (~isempty (varargin))
    for ii=1:2:length(varargin)-1
        if (strcmpi (varargin {ii}, 'value'))
            value = varargin {ii+1};
        elseif (strcmpi (varargin {ii}, 'gradient'))
            gradient = varargin {ii+1};
        elseif (strcmpi (varargin {ii}, 'hessian'))
            hessian = varargin {ii+1};
        else
            error ('sp_tile: unknown option %s', varargin {ii});
        end
    end
end

ndim = size (msh_tile.quad_pnts_in_srf_sp, 1);

knots = space.space_untrimmed.knots;
degree = space.space_untrimmed.degree;

% Nasty thing. We should have two different variables: the dimension of the
% parametric domain (ndim_geo) and the dimension of the domain on which we
% integrate (ndim_quad)
ndim_aux = msh_tile.ndim; msh_tile.ndim = ndim;
sp_on_trimmed = sp_evaluate_element_list (space.space_untrimmed, msh_tile, 'value', false);
msh_tile.ndim = ndim_aux;

sp_tile = struct('ncomp', sp_on_trimmed.ncomp, 'nsh_max', sp_on_trimmed.nsh_max, 'nsh', ... 
            sp_on_trimmed.nsh, 'ndof', space.ndof, 'ndof_dir', sp_on_trimmed.ndof_dir, ...
            'connectivity', space.global_to_active(sp_on_trimmed.connectivity));

% computing the shape functions and their derivatives if needed
if (value || gradient || hessian)
    for idim = 1:ndim
        qn = reshape (msh_tile.quad_pnts_in_srf_sp(idim,:,:), msh_tile.nqn, msh_tile.nel);
        s = findspan(numel(knots{idim})-degree(idim)-2, degree(idim), qn, knots{idim});
        ders_0 = basisfunder (s, degree(idim), qn, knots{idim}, 0);
        ders_value{idim} = permute (reshape(ders_0, [msh_tile.nqn, msh_tile.nel, degree(idim)+1]), [1 3 2]);
    end
    if (value)
        sp_tile.shape_functions = zeros(msh_tile.nqn, sp_on_trimmed.nsh_max, msh_tile.nel);
    end
end

if (gradient || hessian)
    for idim = 1:ndim
        qn = msh_tile.quad_pnts_in_srf_sp(idim,:,:);
        s = findspan(numel(knots{idim})-degree(idim)-2, degree(idim), qn, knots{idim});
        ders_1 = basisfunder (s, degree(idim), qn, knots{idim}, 1);
        ders_1 = ders_1(:, 2, :);

        ders_grad{idim} = permute (reshape(ders_1, [msh_tile.nqn, msh_tile.nel, degree(idim)+1]), [1 3 2]);
    end
    if (gradient)
        sp_tile.shape_function_gradients = zeros(ndim, msh_tile.nqn, sp_on_trimmed.nsh_max, msh_tile.nel);
    end
end

if (hessian) % for now implemented up to the second derivatives
    for idim = 1:ndim
        qn = msh_tile.quad_pnts_in_srf_sp(idim,:,:);
        s = findspan(numel(knots{idim})-degree(idim)-2, degree(idim), qn, knots{idim});
        ders_2 = basisfunder (s, degree(idim), qn, knots{idim}, 2);
        ders_2 = ders_2(:, 3, :);
        ders_hess{idim} = permute (reshape(ders_2, [msh_tile.nqn, msh_tile.nel, degree(idim)+1]), [1 3 2]);
    end  
    sp_tile.shape_function_hessians = zeros(ndim, ndim, msh_tile.nqn, sp_on_trimmed.nsh_max);
end

%%% Computations only for the case ndim = 2
if (value)
  for iel = 1:msh_tile.nel
    for iq = 1:msh_tile.nqn
      sp_tile.shape_functions(iq,:,iel) = kron(ders_value{2}(iq,:,iel), ders_value{1}(iq,:,iel));
    end
  end
end
    
if (gradient)
  for iel = 1:msh_tile.nel
    for iq = 1:msh_tile.nqn
      sp_tile.shape_function_gradients(1, iq, :, iel) = kron(ders_value{2}(iq,:,iel), ders_grad{1}(iq,:,iel)); % computing the tensor product shape functions 
      sp_tile.shape_function_gradients(2, iq, :, iel) = kron(ders_grad{2}(iq,:,iel), ders_value{1}(iq,:,iel)); % computing the tensor product shape functions 
    end
  end
end
    
if (hessian)
  for iel = 1:msh_tile.nel
    for iq = 1:msh_tile.nqn
      sp_tile.shape_function_hessians(1, 1, iq, :, iel) = kron(ders_value{2}(iq,:,iel), ders_hess{1}(iq,:,iel));
      sp_tile.shape_function_hessians(1, 2, iq, :, iel) = kron(ders_grad{2}(iq,:,iel), ders_grad{1}(iq,:,iel));
      sp_tile.shape_function_hessians(2, 1, iq, :, iel) = sp_tile.shape_function_hessians(1, 2, iq, :,iel);
      sp_tile.shape_function_hessians(2, 2, iq, :, iel) = kron(ders_hess{2}(iq,:,iel), ders_value{1}(iq,:,iel));
    end
  end
end

end
