function [r] = cutest_hprod_STRU_AL( x, v, p )
%------------ augmented Lagrangian structured Hessian vector product------%
%
% cutest_hprod_STRU.m is a function to compute the Hessian vector product
% for structured Hessian matrices, when defined by augmented Lagrangians. 
% The augmented Lagrangian is:
%
% L = f(x) + 0.5MU_AL* c(x)^T*c(x),
%
% where f is the CUTEst objective (known Hessian) and c(x)^T*c(x) has
% an unknown Hessian. 
%
%-------------------------------------------------------------------------%
% 09/10/19, J.B., initial version

  r = cutest_hprod(x, v, p);
  
  return;
  
end