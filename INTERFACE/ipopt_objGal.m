function [ gval ] = ipopt_objGal( x, mu )
%------------ unconstrained augmented Lagrangian objective ---------------%
%
% ipopt_objGal.m is a function to compute the Augmented Lagrangian
%
% L = f(x) + 0.5mu c(x)^T*c(x),
%
% where f is the CUTEst objective and c(x^*) = 0.
%
% INPUT:
% x := Current iterate.
% mu:= Penalty paramter.
%
% OUTPUT:
% gval := Gradient value.
%
%-------------------------------------------------------------------------%
% Initial version: 08/23/19, J.B.

c   = cutest_cons(x);
%[c,J]   = cutest_cons(x);

% Gradient evaluation
[~,g]   = cutest_obj(x);

%gval    = g + mu.*J'*c;
gval    = g + mu.*cutest_jtprod(x,c);

end

