function l2_boundary_norm = get_l2_norm_boundary(boundaries,space,msh,gamma_side,problem_data,alpha_vect)
%compute the l2 norm squared  of the
%function $-\alpha\frac{log(\sqrt{x^2+y^2}){2\pi}$. The integration
%boundary has to be specified by the user.

Nbnd = cumsum ([0, boundaries.nsides]);
l2_boundary_norm = 0;
d_gamma = @( x , y , ind ) problem_data.h( x, y, ind ) - problem_data.u_0_exact( x , y , alpha_vect) ; 

for iref = gamma_side
      iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
      d_gamma_squared = @(varargin) ( d_gamma(varargin{:},iref) ) .^ 2;
      l2_boundary_norm = l2_boundary_norm + line_int_d_gamma_mp_defeaturing (space.boundary, msh.boundary, d_gamma_squared, iref_patch_list);
      
end

return