function [ cvals ] = ipopt_cons( x ) % , mu
%-------------------------- ipopt constraints  ---------------------------%
%
% ipopt_cons.m is a function to interface with the CUTEst constraints
%
% c(x).
%
%
% INPUT:
% x := Current iterate.
%
% OUTPUT:
% cvals := Constraint values.
%
%-------------------------------------------------------------------------%
% Initial version: 09/06/19, J.B.

% Constraint evaluations
cvals   = cutest_cons(x);

return;

end

