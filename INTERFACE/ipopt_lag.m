function [ lagval ] = ipopt_lag(x,sigma,lambda)
%-------------------------- IPOPT Lagrangian -----------------------------%
%
% ipopt_lag.m is a function to interface with CUTEst to
% compute the Lagrangian.
%
% L = sigma.* f''(x) + c''(x)*lambda
%
% This function modifies the Lagrangian returned by CUTEst to obtain 
% the form of L.
%
% INPUT:
% x     := Current iterate.
% sigma := IPOPT scaling.
% lambda:= Lagrange multipliers.
%
% OUTPUT:
% lagval := Lagrangian for IPOPT.
%
%-------------------------------------------------------------------------%
% Initial version: 09/06/19, J.B.

% CUTEst objective Hessian
sh      = sparse(cutest_ihess(x,0));

% CUTEst Lagrangian
sl      = sparse(cutest_hess(x,lambda(:)));

lagval  = tril(sl+(sigma-1).*sh);

return;

end

