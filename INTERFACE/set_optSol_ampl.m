% unconstrainted nlp solver
%
%- Nai-Yuan Chiang 2015
function [prob]=set_optSol_ampl(name,prob_in)
prob = prob_in;
[x0,xL,xU,lam0,cL,cU] = spamfunc(name);
% opt solution
        prob.xsol = x0;   
        
end