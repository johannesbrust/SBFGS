%%%%%%% line-search CSBM1 (compact-structured-BFGS Minus Version 1) %%%%%%%
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. Structured objective functions are assumed to be
%
% f(x) = k(x) + u(x), x \in R^n,
%
% where k(x) has a known Hessian, i.e., k''(x) = K(x), and u(x) has an 'unknown'
% Hessian, i.e., u''(x) = U(x). This function extends the original
% line-search method.
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leffeyer, M.Anitescu.
% A report on the compact representations is in DOCS/report_compact_sQN_011819 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 04/05/19
% 04/08/19, J.B., modification to use line-search cvsrch_CSBM1, which
% returns Hessians and vectors bk computed during the line-search.


% REQUIRE
% f0
% g0
% p (step)

function [ret,newIter,newStep,newProb,num_LS,use_zoom,par] = line_search_CSBM1(iter,step,prob,par_in,func,options)

%% Line-search options
c1 			= options.c1;
c2 			= options.c2;

alp_max   	= options.alp_max;
alp_init   	= options.alp_init;

xtol        = options.xtol;
stpmin      = options.stpmin;
stpmax      = options.stpmax;
maxfev      = options.maxfev;

alg			= par_in.alg;
par         = par_in;

n           = length(iter.x);
x           = iter.x;

f0          = prob.obj;
g0          = prob.obj_grad;
p           = step.x;


% line-search info
num_LS = 1;

% check if direction is descent

% if g0'*p >=0
%     aws=g0'*p;
%     warning('Not decsent direction');
% end

% start line search loop
%doLoop  = 1;
%maxflag = 0;

findwolfe   = 0;
zoomFlag    = 1; % ? Use of flag
ret2        = 1;

use_zoom    = 0;
par.skipBFGS = 0;

% Line-search computes, bk, gamma, fWol, gWol
% Unless info==3 (max. evaluations reached)
% it should return points with bk'sk > 0. Here out_CSBM1 contains fields:
%   csbm1_out.bk    = ybar + Kk*sk
%   csbm1_out.ybar  = ybar
%   csbm1_out.Hes   = Hess(x_{k1})
%   csbm1_out.gamma = bk'sk
%   csbm1_out.nhes  = # Hessian evaluations
%   csbm1_out.sk    = x_{k1} - xk
[xWol,fWol,gWol,alp_Wol,info,nfev, stx, sty, out_CSBM1] = cvsrch_CSBM1(@fcn,n,x,f0,g0,p,alp_init,c1,c2,xtol,stpmin,stpmax,maxfev,prob,par,func);

if info==0
    warning('Improper input parameters.');
elseif info==1
    findwolfe=1;
elseif info==2
    warning('Relative width of the interval of uncertainty is too tight.');
elseif info==3
    warning('Number of function evaluation has reached its limit.');
    maxflag=1;
elseif info==4
    warning('The step is at the lower bound stpmin.');
    % 	  findwolfe=1;
elseif info==5
    warning('The step is at the upper bound stpmax.');
    findwolfe=1;
elseif info==6
    warning('Rounding errors prevent further progress.');
elseif info==7
    %  a=[0.0:0.01:0.1]; hold on; for (ii=1:length(a)) px=iter.x+a(ii)*p;   fct(ii) = fcn(px,prob,par,func); end;    plot(a,fct,'-');
    warning('Will try with (much) larger upper bound.');
%     fprintf('stx %g   sty %g\n', stx,sty);

    [xWol,fWol,gWol,alp_Wol,info,nfev,stx,sty, out_CSBM1] = cvsrch_CSBM1(@fcn,n,x,f0,g0,p,alp_init,c1,c2,xtol,stpmin,4000,maxfev,prob,par,func);
    findwolfe=1;
    
    if info==7
        warning('Second line search failed. Function is close to linear.');
        
        par.skipBFGS = 1;
        alp_Wol = 100;
        xWol = x + alp_Wol*p;
        [fWol,gWol] = fcn(xWol,prob,par,func);
        nfev = nfev + 1;
        
    elseif info~=1
        warning('Second line search failed.');
        
        findwolfe=0;
        
    end
end

par.zoomLS = par.zoomLS+nfev;


trialIter    		= iter;
trialIter.x			= xWol;

trialProb			= prob;
trialProb.obj 		= fWol;
trialProb.obj_grad 	= gWol;

trialStep     		= step;
trialStep.size  	= alp_Wol;

%% 04/08/19, J.B., store data from line-search

trialProb.bk        = out_CSBM1.bk;
trialProb.yk        = out_CSBM1.yk;
trialProb.Hes       = out_CSBM1.Hes;
trialProb.gamma     = out_CSBM1.gamma;
trialProb.nhes      = out_CSBM1.nhes;
trialProb.sk        = out_CSBM1.sk;

% find correct step-size for s-bfgs method
if (findwolfe==1)
    
    trialProb  		= evalFunc(trialIter,trialProb,par,func,2);
    
    if strcmp(alg,'s-bfgs')  && zoomFlag==1 && par.skipBFGS == 0 && ...
            trialProb.gamma < 0 % Extra check
        
        [ret2,trialProb,par,couldnotfindtildealpha] = modify_WolfeStep(iter,prob,trialIter,trialStep,trialProb,par,func, stx, sty);
        
        if couldnotfindtildealpha 
            warning('this should not happen');
            ret = 0;
            return;
        end
    end
end

newIter    	= trialIter;
newProb		= trialProb;
newStep     = trialStep;

if findwolfe==1 && ret2==1
    ret = 1;
else
    ret=0;
end

end


function [f,g,gradu, Hkplus] = fcn(x,prob,par,func)
trialIter.x = x;
trial_prob  = evalFunc(trialIter,prob,par,func,2);
f = trial_prob.obj;
g = trial_prob.obj_grad;
gradu = trial_prob.obj_grad_u;
Hkplus = trial_prob.Hes;
end


