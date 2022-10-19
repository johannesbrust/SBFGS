function [ fval ] = ipopt_objal( x, mu ) % , mu
%------------ unconstrained augmented Lagrangian objective ---------------%
%
% ipopt_objal.m is a function to compute the Augmented Lagrangian
%
% L = f(x) + 0.5mu c(x)^T*c(x),
%
% where f is the CUTEst objective and c(x^*) = 0. This is intended to
% use with IPOPT.
%
% INPUT:
% x := Current iterate.
% mu:= Penalty paramter.
%
% OUTPUT:
% fval := Function value.
%
%-------------------------------------------------------------------------%
% Initial version: 08/23/19, J.B.

% Function evaluation
f       = cutest_obj(x);
c       = cutest_cons(x);

fval    = f + 0.5*mu*(c'*c);

end

