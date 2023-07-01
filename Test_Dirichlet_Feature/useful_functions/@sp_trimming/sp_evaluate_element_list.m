% SP_EVALUATE_ELEMENT_LIST: compute the basis functions in a given list of elements.
%
%     sp = sp_evaluate_element_list (space, msh_elems, 'option1', value1, ...)
%
% INPUTS:
%     
%    space:     object defining the space of discrete functions (see sp_trimming)
%    msh_elems: msh structure containing the information of quadrature or
%               visualization points, for a given list of elements 
%               (see msh_cartesian/msh_evaluate_element_list)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+----------------------------------
%            value      |      true       |  compute shape_functions
%            gradient   |      false      |  compute shape_function_gradients
%            hessian    |      false      |  compute shape_function_hessians
%            laplacian  |      false      |  compute shape_function_laplacians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                     DESCRIPTION
%    ncomp           (scalar)                                   number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                                   total number of degrees of freedom
%    ndof_dir        (1 x ndim vector)                          degrees of freedom along each direction
%    nsh_max         (scalar)                                   maximum number of shape functions per element
%    nsh             (1 x msh_elems.nel vector)                 actual number of shape functions per each element
%    connectivity    (nsh_max x msh_elems.nel vector)           indices of basis functions that do not vanish in each element
%    shape_functions (msh_elems.nqn x nsh_max x msh_elems.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (rdim x msh_elems.nqn x nsh_max x msh_elems.nel)        basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%        rdim x rdim x msh_elems.nqn x nsh_max x msh_elems.nel) basis function hessians evaluated at each quadrature node in each element
%    shape_function_laplacians 
%       (msh_elems.nqn x nsh_max x msh_elems.nel)               basis functions laplacians evaluated at each quadrature node in each element
%
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

function sp = sp_evaluate_element_list (space, msh, varargin)

sp = sp_evaluate_element_list (space.space_untrimmed, msh, varargin{:});
sp.connectivity = space.global_to_active (sp.connectivity);
sp.ndof = space.ndof;

end
