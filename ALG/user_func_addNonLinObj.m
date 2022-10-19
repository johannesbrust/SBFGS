%%% user-provided functions
function func = user_func_addNonLinObj(rawfunc)

%%%%%% functions
 func.rawfunc    	= rawfunc;
 func.getObj    	= @getObj;
 func.getGrad    	= @getGrad;
 func.getKnowGrad   = @getKnowGrad;
 func.getUnknowGrad = @getUnknowGrad;
 func.getHess   	= @getHess;


%%%%%% get objective function
function  f = getObj(iters,prob,par,rawfunc)
f=0;
f = rawfunc.getRawObj(iters.x);
% add extra term
for k = 1:par.num_unknown
    f = f + log((iters.x(k)-prob.xsol(k))^2+1);
end

%%%%%% get gradient 
function [g,gk,gu] = getGrad(iters,prob,par,rawfunc)
gk = getKnowGrad(iters,prob,par,rawfunc);
gu = zeros(length(iters.x),1);
for k = 1:par.num_unknown
    gu(k) = 2*(iters.x(k)-prob.xsol(k))/((iters.x(k)-prob.xsol(k))^2+1);
end
g = gk+gu;

%%%%%% get gradient of known part
function gk = getKnowGrad(iters,prob,par,rawfunc)
gk = rawfunc.getRawGrad(iters.x);

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iters,prob,par,rawfunc)
gu = zeros(length(iters.x),1);
for k = 1:par.num_unknown
    gu(k) = 2*(iters.x(k)-prob.xsol(k))/((iters.x(k)-prob.xsol(k))^2+1);
end

% get Hessian ( known part only )
function H = getHess(iters,prob,par,rawfunc)
H = rawfunc.getRawHess(iters.x);
if strcmp(par.alg,'exact')
    for k = 1:par.num_unknown
        H(k,k) = H(k,k) + (2-2*(iters.x(k)-prob.xsol(k))^2)/((iters.x(k)-prob.xsol(k))^2+1)^2;
    end  
end












