function [r] = QUAD_hprod_STRU(s, Q1)
%------------ QUAD_hprod_STRU structured Hessian product -----------------%
%
% The Hessian product with a structured quadratic Hessian:
%
% r = Q1*s,
%
%-------------------------------------------------------------------------%
% 10/23/19, J.B.

  r = Q1*s;
  
  return;
  
end