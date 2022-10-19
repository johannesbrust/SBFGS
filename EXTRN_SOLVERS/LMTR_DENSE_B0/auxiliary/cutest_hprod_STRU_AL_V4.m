function [r] = cutest_hprod_STRU_AL_V4( x, v, p )
%------------ augmented Lagrangian structured Hessian vector product V4 --%
%
% cutest_hprod_STRU_AL_V4.m is a function to compute the Hessian vector product
% for structured Hessian matrices, when defined by an modified augmented Lagrangian. 
% The augmented Lagrangian is:
%
% L = f(x) + 0.5MU_AL* c(x)^T*c(x),
%
% where f is the CUTEst objective (known Hessian) and c(x)^T*c(x) has
% an unknown Hessian. 
%
%-------------------------------------------------------------------------%
% 09/10/19, J.B., initial version

  global MU_AL;

  r     = cutest_hprod(x, v, p);  
  %c     = cutest_cons(x);  
  
  r1    = MU_AL.*cutest_jprod(x,p);  
  r     = r + cutest_jtprod(x,r1);
  
  return;
  
end