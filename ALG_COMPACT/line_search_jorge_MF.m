%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% line search method for wolfe condition
%
%

function [ret,newIter,newStep,newProb,num_LS,use_zoom,par] = line_search_jorge_MF(iter,step,prob,par_in,func)

alg			= par_in.alg;
par         = par_in;
bfgsAlg		= par.sbfgsAlg;


alp_max   	= 2.0;
alp_init   	= 1.0;

f0		= prob.obj;
g0		= prob.obj_grad;
p		= step.x;

% alg parameters
maxLSstep		= par.maxLSstep;
c1 				= par.c1;
c2 				= par.c2;

% line-search info
num_LS = 1;

% check if direction is descent
if g0'*p >=0
    aws=g0'*p;
    warning('Not decsent direction');
end

% start line search loop
doLoop  = 1;
maxflag = 0;
findwolfe=0;
zoomFlag = 1;
ret2 = 1;

use_zoom = 0;
par.skipBFGS = 0;

n = length(iter.x);
x = iter.x;
xtol = 1e-6;
stpmin = 0;
stpmax = alp_max;
maxfev	= 3000;
nfev    = 0;
[xWol,fWol,gWol,alp_Wol,info,nfev, stx, sty] = cvsrch_MF(@fcn,n,x,f0,g0,p,alp_init,c1,c2,xtol,stpmin,stpmax,maxfev,prob,par,func);

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
    [xWol,fWol,gWol,alp_Wol,info,nfev,stx,sty] = cvsrch_MF(@fcn,n,x,f0,g0,p,alp_init,c1,c2,xtol,stpmin,4000,maxfev,prob,par,func);
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

% find correct step-size for s-bfgs method
if (findwolfe==1)
    trialProb  		= evalFunc(trialIter,trialProb,par,func,2);
    if strcmp(alg,'s-bfgs')  && zoomFlag==1 && par.skipBFGS == 0 
        [ret2,trialProb,par,couldnotfindtildealpha] = modify_WolfeStep_MF(iter,prob,trialIter,trialStep,trialProb,par,func, stx, sty);
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


