% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    object defining the space of discrete functions (see sp_trimming)
%  msh:   object defining the domain partition and the quadrature rule (see msh_trimming)
%  h:     function handle to compute the Dirichlet condition
%  sides: labels of boundaries on which a Dirichlet condition is imposed.
%           The boundary sides cannot be on the trimming curves.
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: numbering of the corresponding basis functions, in the trimmed space
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

function [u, dofs] = sp_drchlt_l2_proj (space, msh, h, sides)

  param_sides = [];
  for side = sides
      bnd_id = msh_trimming_boundary_find (msh.reparam, side);
      param_side = msh.reparam.boundaries(bnd_id).parametric_side;
      if param_side > 0
          param_sides = [param_sides, param_side];
      else
          error('Strong Dirichlet boundary conditions must be applied to non trimmed boundaries!');
      end
  end

  [u, dofs] = sp_drchlt_l2_proj (space.space_untrimmed, msh.msh_cart, h, param_sides);
  
  [~, position, dofs] = intersect (dofs, space.active_dofs);
  u = u(position);
  
end
