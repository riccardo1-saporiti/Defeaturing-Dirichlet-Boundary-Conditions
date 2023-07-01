% SP_EVALUATE_TILES: create a struct of a scalar space defined over tiles coming from the reparametrization of a trimmed element,
%                                       where the evaluation is performed in the physical domain
% 
%   sp_tile = sp_evaluate_tiles (space, msh_tiles, varargin)
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
%            laplacian  |      false      |  compute shape_function_laplacians%
% OUTPUT:
%
%    sp_tile: struct representing the discrete function space, with the fields described in sp_evaluate_element_list
%              (see the article for a detailed description)
%  The functions will be numbered as in the tensor-product space. This
%    numbering must be changed afterwards.
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

function  sp_tile = sp_evaluate_tiles (space, msh_tiles, varargin)
value = true;
gradient = false;
hessian = false;
laplacian = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_tiles: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'laplacian'))
      laplacian = varargin {ii+1};
    else
      error ('sp_evaluate_tiles: unknown option %s', varargin {ii});
    end
  end
end

hessian_param = hessian || laplacian;
grad_param = gradient || hessian_param;
value_param = value || grad_param;

sp_tile = sp_evaluate_tiles_param (space, msh_tiles, 'value', value_param, 'gradient', grad_param, 'hessian', hessian_param);

switch (lower (space.space_untrimmed.transform))
  case {'grad-preserving'}
    sp_tile = sp_tile_scalar_grad_preserving_transform(sp_tile, msh_tiles, value, gradient, hessian, laplacian);
  case {'integral-preserving'}
    sp_tile = sp_tile_scalar_integral_preserving_transform(sp_tile, msh_tiles, value);
    if (gradient || hessian || laplacian)
      warning ('Derivatives not implemented for integral-preserving transformation')
    end
end

end