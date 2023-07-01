function errh1s = trimming_errh1s_negative (u_0, msh_cart, space, reparam, u_ex)

A_gradu_gradv = op_gradu_gradv_trimming(space, space, msh_cart, reparam);
errh1s = sqrt((u_0-u_ex)'*A_gradu_gradv*(u_0-u_ex));

end