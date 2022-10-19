%%% user-provided functions
function func = user_func_setRatio(rawfunc)

%%%%%% functions
 func.rawfunc    	= rawfunc;
 func.getObj    	= @getObj;
 func.getGrad    	= @getGrad;
 func.getKnowGrad   = @getKnowGrad;
 func.getUnknowGrad = @getUnknowGrad;
 func.getHess   	= @getHess;


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

% get Hessian ( known part only )
function H = getHess(iter,prob,par,rawfunc)
H = par.exactRatio * rawfunc.getRawHess(iter.x);












