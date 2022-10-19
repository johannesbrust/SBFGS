function [ gval ] = ipopt_objG( x )
%-------------------------- unconstrained objective ----------------------%
%
% ipopt_objG.m is a function to interface with the CUTEst objective
%
% f(x).
%
%
% INPUT:
% x := Current iterate.
%
% OUTPUT:
% gval := Gradient value.
%
%-------------------------------------------------------------------------%
% Initial version: 09/06/19, J.B.

% Gradient evaluation
[~,gval]   = cutest_obj(x);

return;

end

