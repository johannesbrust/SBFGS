%%%%%%%%%%%%%%% CS-BFGS (User Function Set Ratio Sparse) %%%%%%%%%%%%%%%%%
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. The compact forms are
%
% B = B0 - Psi M Psi',
%
% where typically B0 represents an initial matrix, and
% Psi (n x 2m), M (2m x 2m) are rectangular and small square matrices,
% correspondingly. 
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leffeyer, M.Anitescu.
% An initial report on the compact representations is in DOCS/report_compact_sQN_011819 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 04/09/19
function func = user_func_setRatio_sp(rawfunc)

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
% Additional call to sparse
function H = getHess(iter,prob,par,rawfunc)
H = par.exactRatio * sparse(rawfunc.getRawHess(iter.x));












