%%%%%%%%%%%%% csBFGS (compact-structured-BFGS Methods) %%%%%%%%%%%%%%%%%%%%
% read_cutest_prob_cons is a function that reads constrained CUTEst
% problems.
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. The compact forms are
%
% B = B0 - Psi M Psi',
%
% where typically B0 = gamma.I (n x n) multiple of identity initial matrix,
% Psi (n x 2m), M (2m x 2m) are small low rank updates. 
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leyffer, M.Anitescu.
%
% This script modifies the original 'read_ampl_prob' from Nai-Yuan 2015.
% Specifically, this function uses CUTEst instead of AMPL.
%
% This function assumes that it is called from within the 'INTERFACE/'
% folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 09/09/19

function [ret, prob, par, rawfunc]=read_cutest_prob_cons(name,prob_in,par_in)

par  = par_in;
prob = prob_in;

% Build CUTEst problem
%probname    = 'BOX';
%cmdcutest   = ['cutest2matlab_osx $MASTSIF/' name];

% Modified build script to access sparse arrays.
cmdcutest   = ['cutest2matlab_osx $MASTSIF/' name]; 
unix(cmdcutest);

% Problem data
cu_prob     = cutest_setup();

% cu_grad     = @cutest_grad;
% cu_obj      = @cutest_obj;

x0          = cu_prob.x;

% initial call to ampl model
%[x0,xL,xU,lam0,cL,cU] = spamfunc(name);


ret = 1;
% check if input is unconstrainted 
% if length(lam0)>0
%   ret =0; 
%   warning('require unconstrainted prob');
% end
% 
% for i = 1:length(x0)
%   if xL(i)>=-1e+20
%     ret =0;
% 	warning('require unconstrainted prob');
%   elseif xU(i)<=1e+20
%     ret =0;
% 	warning('require unconstrainted prob');
%   end
% end

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
% rawfunc.getRawObj = @getRawObj;
% rawfunc.getRawGrad = @getRawGrad;
% rawfunc.getRawHess = @getRawHess;
    
rawfunc.getRawObj   = @cutest_obj;
rawfunc.getRawGrad  = @cutest_grad;
rawfunc.getRawHess  = @cutest_hess;

% Constraint data
prob.m = cu_prob.m;
prob.v = cu_prob.v;
prob.nnzh = cu_prob.nnzh;
prob.nnzj = cu_prob.nnzj;

%------------------ user-provided functions

% get objective function
% function  f = getRawObj(x);
% 
% [f,cons] = spamfunc(x,0);
% 
%   
% % get gradient of objective
% function g = getRawGrad(x)
% 
%   [f,c] = spamfunc(x,0);
%   [g,J] = spamfunc(x,1);
% 
%        
% % get Hessian of objective (wrt x and lam) - assembled
% function H = getRawHess(x)
% 
%   [f,c] = spamfunc(x,0);
%   [g,J] = spamfunc(x,1); 
%   lam_temp = zeros(0,1);
%       H = spamfunc(lam_temp);     
