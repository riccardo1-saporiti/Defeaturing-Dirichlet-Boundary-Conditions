% SP_H1_ERROR: Evaluate the error in L^2 norm in a trimmed domain.
%
%   errl2 = sp_l2_error (space, msh, u, uex, trimmed_elems, non_trimmed_elem_ids);
%
% INPUT:
%
%     space:    object defining the space of discrete functions (see sp_trimming)
%     msh:      object defining the untrimmed domain partition and the quadrature rule (see msh_trimming)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     reparam: a structure containing all the info coming from the reparametrization tool 
%
% OUTPUT:
%
%     errl2:  error in L^2 norm
%
% Copyright (C) 2020 Ondine Chanon, Luca Coradello, Riccardo Puppi
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function errl2 = sp_l2_error (space, msh_trimmed, u, uex)

u_glob = zeros (space.space_untrimmed.ndof, 1);
u_glob(space.active_dofs) = u;
% untrimmed part
errl2_untr = 0; 
if not(isempty(msh_trimmed.reparam.non_trimmed_elem_ids))
    msh_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids);
    sp = sp_evaluate_element_list(space.space_untrimmed, msh_elems, 'value', true);
    
    errl2_untr = sp_l2_error (sp, msh_elems, u_glob, uex);
end

% trimmed part
msh_tiles = msh_evaluate_tiles (msh_trimmed);
sp_tiles = sp_evaluate_tiles (space, msh_tiles, 'value', true);

errl2_tr = sp_l2_error (sp_tiles, msh_tiles, u, uex);

errl2 = sqrt (errl2_tr^2 + errl2_untr^2);

end
