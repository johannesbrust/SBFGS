% unconstrainted nlp solver
%
%- Nai-Yuan Chiang 2015
function [ret, prob, par, rawfunc]=read_ampl_prob(name,prob_in,par_in)

par  = par_in;
prob = prob_in;

% initial call to ampl model
[x0,xL,xU,lam0,cL,cU] = spamfunc(name);
ret = 1;
% check if input is unconstrainted 
if length(lam0)>0
  ret =0; 
  warning('require unconstrainted prob');
end

for i = 1:length(x0)
  if xL(i)>=-1e+20
    ret =0;
	warning('require unconstrainted prob');
  elseif xU(i)<=1e+20
    ret =0;
	warning('require unconstrainted prob');
  end
end

% problem dimensions
       prob.n = length(x0);      % number of variables

% adjust par.addUnknown
if(par.num_unknown > prob.n)
    par.num_unknown = 1;
end
       
% initial point 
if par.initX0 == 1
        prob.x0 = ones(prob.n,1);
else
        prob.x0 = x0;
end

    
% functions
    rawfunc.getRawObj = @getRawObj;
rawfunc.getRawGrad = @getRawGrad;
   rawfunc.getRawHess = @getRawHess;
    
%------------------ user-provided functions

% get objective function
function  f = getRawObj(x);

[f,cons] = spamfunc(x,0);

  
% get gradient of objective
function g = getRawGrad(x)

  [f,c] = spamfunc(x,0);
  [g,J] = spamfunc(x,1);

       
% get Hessian of objective (wrt x and lam) - assembled
function H = getRawHess(x)

  [f,c] = spamfunc(x,0);
  [g,J] = spamfunc(x,1); 
  lam_temp = zeros(0,1);
      H = spamfunc(lam_temp);     
