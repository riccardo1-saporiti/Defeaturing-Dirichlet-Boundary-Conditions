%Compute the line integral I = \int_{\gamma} f ds, exploiting the tensor
%product structure 
%   int = line_int_d_gamma_tp_defeaturing (space, msh, coeff);
%
% INPUT:
%     
%   spv:   object representing the function space (see sp_scalar)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   coeff: function handle to compute the source function
%
% OUTPUT:
%
%   int: computed integral 
% 
% Copyright (C) 2011, 2017 Rafael Vazquez, edited by Riccardo Saporiti
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

function int = line_int_d_gamma_tp_defeaturing (space, msh, coeff)

  for idim = 1:msh.ndim
    size1 = size (space.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('The discrete space is not associated to the mesh')
    end
  end
  
  int = 0 ; 
 
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);

    for idim = 1:msh.rdim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    int = int + line_int_d_gamma_defeaturing (sp_col, msh_col, coeff (x{:}));
  end

end
