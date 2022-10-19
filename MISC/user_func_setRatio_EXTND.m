function func = user_func_setRatio_EXTND(rawfunc)
% user_func_setRatio_EXTND: Extension of setRatio function
%
% INPUT:
% rawfunc   := CUTEst fuctionality
%
% OUTPUT:
% func      := Function struct for algorithms
%-------------------------------------------------------------------------%
% Initial version: 11/19/19, J.B.

%%%%%% functions
 func.rawfunc    	= rawfunc;
 func.getObj    	= @getObj;
 func.getGrad    	= @getGrad;
 func.getKnowGrad   = @getKnowGrad;
 func.getUnknowGrad = @getUnknowGrad;
 func.getHess   	= @getHess;
 func.getHessFull   = @(x,par)(getHessFull(x,par,rawfunc));
 func.getHessProd   = @(x,s,par)(getHessProd(x,s,par,rawfunc)); 


%%%%%% get objective function
function  f = getObj(iter,prob,par,rawfunc)
f = 0;
f = rawfunc.getRawObj(iter.x);

%%%%%% get gradient 
function [g,gk,gu] = getGrad(iter,prob,par,rawfunc)
g  = rawfunc.getRawGrad(iter.x);
gk = par.exactRatio * g;
gu = g-gk;

%%%%%% get gradient of known part
function gk = getKnowGrad(iter,prob,par,rawfunc)
gk = par.exactRatio * rawfunc.getRawGrad(iter.x);

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iter,prob,par,rawfunc)
gu = (1-par.exactRatio) * rawfunc.getRawGrad(iter.x);

%% Extended functions
% get Hessian ( known part only )
function H = getHess(iter,prob,par,rawfunc)
H = 0;

    if strcmp(par.alg,'s-bfgs') == true &&...
            par.withHess == true
        
    H = par.exactRatio * rawfunc.getRawHess(iter.x);
        
    end

% getHessFull ( known part only )
function H = getHessFull(x,par,rawfunc)
H = par.exactRatio * rawfunc.getRawHess(x);


% getHessProd ( known part only )
function Hs = getHessProd(x,s,par,rawfunc)
 Hs = par.exactRatio * (rawfunc.getRawHess(x)*s);

    

    












