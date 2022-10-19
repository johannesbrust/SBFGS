%%% user-provided functions
function func = user_func_addQuadObj(rawfunc)

%%%%%% functions
 func.rawfunc    	= rawfunc;
 func.getObj    	= @getObj;
 func.getGrad    	= @getGrad;
 func.getKnowGrad   = @getKnowGrad;
 func.getUnknowGrad = @getUnknowGrad;
 func.getHess   	= @getHess;


%%%%%% get objective function
function  f = getObj(iter,prob,par)
f=0;
f = rawfunc.getRawObj(iter.x);
% add extra term
for k = 1:par.num_unknown
    f = f + (x[k]-prob.x0[k])^2;
end

%%%%%% get gradient 
function [g,gk,gu] = getGrad(iters,prob,par)
gk = rawfunc.getRawGrad(iter.x);
% add extra term
gu = zeros(length(iter.x));
for k = 1:par.num_unknown
    gu[k] = 2*(x[k]-prob.x0[k]);
end
g = gk+gu;

%%%%%% get gradient of known part
function gk = getKnowGrad(iters,prob,par)
gk = rawfunc.getRawGrad(iter.x);

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iters,prob,par)
gu = zeros(length(iter.x));
for k = 1:par.num_unknown
    gu[k] = 2*(x[k]-prob.x0[k]);
end

% get Hessian ( known part only )
function H = getHess(iters,prob,par)
H = rawfunc.getRawHess(iter.x);












