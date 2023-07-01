% OP_F_V_MP: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i), in a multipatch domain.
%
%   rhs = op_f_v_mp (spv, msh, coeff, [patches]);
%
% INPUT:
%     
%   spv:     object representing the function space (see sp_multipatch)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   coeff:   function handle to compute the source function
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   rhs: assembled right-hand side
% 
% Copyright (C) 2015 Rafael Vazquez
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

function int = line_int_d_gamma_mp_defeaturing (space, msh, coeff, patch_list)

  if (nargin < 4)
    patch_list = 1:msh.npatch;
  end

  if (space.npatch ~= msh.npatch)
    error ('op_f_v_mp: the number of patches does not coincide')
  end
  
  int = 0;
  for iptc = patch_list
    int_loc = line_int_d_gamma_tp_defeaturing (space.sp_patch{iptc}, msh.msh_patch{iptc}, coeff);
    
    if (~isempty (space.dofs_ornt))
      int_loc = space.dofs_ornt{iptc}(:) .* int_loc(:);
    end
    int = int + int_loc;
  end

end