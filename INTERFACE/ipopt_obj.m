function [ fval ] = ipopt_obj( x ) % , mu
%-------------------------- unconstrained objective ----------------------%
%
% ipopt_obj.m is a function to interface with the CUTEst objective
%
% f(x).
%
%
% INPUT:
% x := Current iterate.
%
% OUTPUT:
% fval := Objective value.
%
%-------------------------------------------------------------------------%
% Initial version: 09/06/19, J.B.

% function evaluation
fval   = cutest_obj(x);

return;

end

