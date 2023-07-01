% SP_H1_ERROR: Evaluate the error in H^1 norm in a trimmed domain.
%
%   [errh1, errl2, errh1s, errh1_elem_untrimmed, errh1_elem_trimmed] = sp_h1_error (space, msh, u, uex, graduex, trimmed_elems, non_trimmed_elem_ids);
%
% INPUT:
%
%     space:    object defining the space of discrete functions (see sp_trimming)
%     msh:      object defining the untrimmed domain partition and the quadrature rule (see msh_trimming)
%     u:        vector of dof weights
%     uex:      function handle to evaluate the exact solution
%     graduex:  function handle to evaluate the gradient of the exact solution
%     reparam: a structure containing all the info coming from the reparametrization tool 
%
% OUTPUT:
%
%     errh1:  error in H^1 norm
%     errl2:  error in L^2 norm
%     errh1s: error in H^1 seminorm
%     errh1_elem_untrimmed:  error in H^1 norm, for each single non-trimmed element
%     errh1_elem_trimmed:  error in H^1 norm, for each single trimmed element
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

function [errh1, errl2, errh1s, errh1_elem_untr, errh1_elem_tr] = sp_h1_error (space, msh_trimmed, u, uex, graduex)

u_glob = zeros (space.space_untrimmed.ndof, 1);
u_glob(space.active_dofs) = u;
% untrimmed part
errl2_untr = 0; errh1s_untr = 0; errh1_elem_untr = [];
if not(isempty(msh_trimmed.reparam.non_trimmed_elem_ids))
    msh_elems = msh_evaluate_element_list(msh_trimmed.msh_cart, msh_trimmed.reparam.non_trimmed_elem_ids);
    sp = sp_evaluate_element_list(space.space_untrimmed, msh_elems, 'gradient', true);
    
    [~, errl2_untr, errh1s_untr, errh1_elem_untr] = sp_h1_error (sp, msh_elems, u_glob, uex, graduex);
end

% trimmed part
msh_tiles = msh_evaluate_tiles (msh_trimmed);
sp_tiles = sp_evaluate_tiles (space, msh_tiles, 'value', true, 'gradient', true);

[errh1_tr, errl2_tr, errh1s_tr, errh1_tiles, errl2_tiles, errh1s_tiles] = sp_h1_error (sp_tiles, msh_tiles, u, uex, graduex);

if (nargout == 5)
  errh1_elem_tr = zeros (1, msh_trimmed.reparam.nb_trimmed_elems);
  element_inds = msh_tiles.trimmed_element_index;
  for iel = 1:msh_trimmed.reparam.nb_trimmed_elems
    errh1_elem_tr(iel) = sum (errh1_tiles(element_inds == iel).^2);
  end
  errh1_elem_tr = sqrt (errh1_elem_tr);
end

errl2 = sqrt (errl2_tr^2 + errl2_untr^2);
errh1s = sqrt (errh1s_tr^2 + errh1s_untr^2);
errh1 = sqrt (errl2^2 + errh1s^2);


end
