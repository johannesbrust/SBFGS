function [ jacval ] = ipopt_jac( x )
%-------------------------- IPOPT Jacobian -------------------------------%
%
% ipopt_jac.m is a function to interface with the CUTEst Jacobian
%
% (c(x))'.
%
%
% INPUT:
% x := Current iterate.
%
% OUTPUT:
% jacval := Jacobian value.
%
%-------------------------------------------------------------------------%
% Initial version: 09/06/19, J.B.

% Constraint evaluation
[~,jacvalf]   = cutest_cons(x);

jacval = sparse(jacvalf);

return;

end

