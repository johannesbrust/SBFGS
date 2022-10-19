function [r] = LIBSVM_hprod_STRU(w)
%------------ LIBSVM structured Hessian product --------------------------%
%
% The Hessian product with a regularized logistic regression Hessian:
%
% r = LAM.w,
%
% where LAM is global regularization parameter.
%-------------------------------------------------------------------------%
% 10/04/19, J.B.

  global LAM;

  r = LAM.*w;
  
  return;
  
end