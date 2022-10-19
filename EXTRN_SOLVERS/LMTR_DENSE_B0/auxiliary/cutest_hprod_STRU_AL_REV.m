function [r] = cutest_hprod_STRU_AL_REV( x, v, p )
%------------ augmented Lagrangian structured Hessian vector product------%
%
% cutest_hprod_STRU_AL_REV.m is a function to compute the Hessian vector product
% for structured Hessian matrices, when defined by augmented Lagrangians. 
% The augmented Lagrangian is:
%
% L = f(x) + 0.5MU_AL* c(x)^T*c(x),
%
% with Hessian:
%
% L'' = f''(x) + MU_AL*(sum(c_i''(x)*c_i) + c'(x)(c'(x))^T).
%
% where f is the CUTEst objective (unknown Hessian) and c(x)^T*c(x) has
% an known Hessian. This function reverses, which Hessian is known or
% unknown and computes the product only with the known part.
%
% v is assumed a vector of zeros.
%-------------------------------------------------------------------------%
% 09/13/19, J.B., initial version

    global MU_AL;

    c  = cutest_cons(x);
    jp = cutest_jprod(x,p);

    % CUTEst Lagrangian product: (f'' + sum(c_i''*c_i))*p
    r  = cutest_hprod(x,c,p);
    
    % Product with only known AL Hessian .
    r  = MU_AL.*(r - cutest_hprod(x,v,p) + cutest_jtprod(x,jp));
  
    return;
  
end