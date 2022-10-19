function func = user_func_AL(rawfunc)
%------------ unconstrained augmented Lagrangian objective ---------------%
%
% user_func_AL.m is a function to compute the Augmented Lagrangian
%
% L = f(x) + 0.5mu c(x)^T*c(x),
%
% where f is the CUTEst objective (known Hessian) and c(x)^T*c(x) has
% an unknown Hessian. The inputs/outputs mimic other user functions
% for the structured BFGS solvers.
%
% INPUT:
% rawfunc := Struct with CUTEst functions.
%
% OUTPUT:
% func := Struct for function evaluations.
%
%-------------------------------------------------------------------------%
% Initial version: 08/23/19, J.B.

func.rawfunc    	= rawfunc;
func.getObj         = @getObj;
func.getGrad    	= @getGrad;
func.getKnowGrad    = @getKnowGrad;
func.getUnknowGrad  = @getUnknowGrad;
func.getHess        = @getHess;

end

%%%%%% get objective function
function  f = getObj(iter,prob,par,rawfunc)

    ff  = cutest_obj(iter.x);
    c   = cutest_cons(iter.x);
    
    f   = ff + 0.5*par.mu*(c'*c);
end

%%%%%% get gradient 
function [g,gk,gu] = getGrad(iter,prob,par,rawfunc)

    [~,gk]      = cutest_obj(iter.x);    
    c           = cutest_cons(iter.x);
    %[c,J]   = cutest_cons(iter.x);
    
    gu  = (par.mu).*cutest_jtprod(iter.x,c);
    %gu      = (par.mu).*J'*c;
    
    g   = gk + gu; 
    
end    
    
%%%%%% get gradient of known part
function gk = getKnowGrad(iter,prob,par,rawfunc)
    gk = cutest_grad(iter.x);
end

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iter,prob,par,rawfunc)

     c   = cutest_cons(iter.x);
     gu = (par.mu).*cutest_jtprod(iter.x,c);
     
     %[c,J]   = cutest_cons(iter.x);
    
    %gu  = (par.mu).*cutest_jtprod(iter.x,c);
    %gu      = (par.mu).*J'*c;
     
end

% get Hessian ( known part only )
function H = getHess(iter,prob,par,rawfunc)
    H = cutest_ihess(iter.x,0);
end
