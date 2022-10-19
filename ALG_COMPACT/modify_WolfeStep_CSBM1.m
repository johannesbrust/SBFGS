%%%%%%%%%%%%%%%%%%%%%% CSBM1 (modify Wolfe Step) %%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. Structured objective functions are assumed to be
%
% f(x) = k(x) + u(x), x \in R^n,
%
% where k(x) has a known Hessian, i.e., k''(x) = K(x), and u(x) has an 'unknown'
% Hessian, i.e., u''(x) = U(x). 
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leffeyer, M.Anitescu.
% A report on the compact representations is in DOCS/report_compact_sQN_011819 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 04/08/19
% 

function [ret,wolfeProb,par, couldnotfindtildealpha] = modify_WolfeStep_CSBM1(currIter,currProb,wolfe_iter,wolfe_step,wolfe_prob,par_in,func, stx,sty)
 couldnotfindtildealpha = 0;
 % alg parameters
 par            = par_in;

 FR				= par.modifyFR;
 maxLSstep		= par.maxLSstep;
 bfgsAlg		= par.sbfgsAlg;
 Pmat           = currProb.Pmat;

 % shortcut
 wolfeIter		= wolfe_iter;
 wolfeStep		= wolfe_step;
 wolfeProb		= wolfe_prob;

 % current BFGS approximation and exact Hessian
 A_k = currProb.Bk;
  
 Hes_k			= currProb.Hes;
 Hes_kplus1		= wolfeProb.Hes;


 % start line search loop  
 alp_max = wolfeStep.size;
 doLoop = 1;
 ret=0;
 updateS = 0;
 num_LS = 1; 
 
 % output 
 sbfgsIter  	= currIter;
 sbfgsStep  	= wolfeStep;
 sbfgsProb		= currProb;

 % try if wolfe step satisfies gamma > 0
 grad_sbfgs		= sbfgsProb.obj_grad_u;
 grad_wolfe		= wolfeProb.obj_grad_u;
 
 if bfgsAlg==1
   	Hes_use =  Hes_k;
 elseif bfgsAlg==2
    Hes_use =  Hes_kplus1;
 else
   	error('bfgsAlg should be 1 or 2');
 end
 
 sk = wolfeIter.x - sbfgsIter.x;
 yk = grad_wolfe  - grad_sbfgs;

 if strcmp(par.addUnknown,'addNonLinObj_Proj')
    sk = Pmat*sk;
    Hes_use = Pmat*Hes_use*Pmat';
    Hes_kplus1 = Pmat*Hes_kplus1*Pmat';
 end 
 
 %ak = (A_k + Hes_use)*sk;
 bk = yk + Hes_kplus1*sk;

 %alpha  = sk'*ak;
 gamma  = sk'*bk;
 gamma_ori =gamma;
 
 % if we can do regularizaion later, we do not care about gamma 
 if gamma>0 || bfgsAlg==2 || (bfgsAlg==1 && strcmp(par.addUnknown,'addNonLinObj_Proj') )
%|| (bfgsAlg==1 && strcmp(par.addUnknown,'addNonLinObj_Proj') )   <---  if we do this, we can solve one more prob
   updateS = 1;
	doLoop = 0;
%  elseif gamma < 0
%      warning(['gamma < 0: ' num2str(gamma_ori)]);
 end
 
 if strcmp(par.addUnknown,'addNonLinObj_Proj') && ( norm(yk,inf) < 1e-16)
	par.stopBFGS=1;
 else 
	if strcmp(par.addUnknown,'addNonLinObj_Proj') && (par.stopBFGS==1)
	  warning('reuse BFGS approximation matrix');	
	end
	par.stopBFGS=0;
 end 

 % line-search info

while doLoop == 1 
   	sbfgsStep.size		= alp_max*FR(num_LS);   
   	%sbfgsIter   	  	= update_iterates(sbfgsIter,sbfgsStep,sbfgsProb,par);
    sbfgsIter   	  	= update_iterates(currIter,sbfgsStep,sbfgsProb,par);

    % compute derivatives
    sbfgsProb  			= evalFunc(sbfgsIter,sbfgsProb,par,func,1);
	grad_sbfgs			= sbfgsProb.obj_grad_u;

	sk = wolfeIter.x - sbfgsIter.x;
    yk = grad_wolfe  - grad_sbfgs;
    
    if strcmp(par.addUnknown,'addNonLinObj_Proj')
         sk = Pmat*sk;
    end 

	ak = A_k*sk + Hes_use*sk;
	bk = yk + Hes_kplus1*sk;

	alpha  = sk'*ak;
	gamma  = sk'*bk;

    %transpose(wolfe_prob.obj_grad-sbfgsProb.obj_grad)
    %yk'
    %transpose(Hes_kplus1*sk)
    %transpose(wolfe_prob.obj_grad-sbfgsProb.obj_grad)*sk
	if gamma>0
	  updateS = 1;
    elseif gamma<-1e-15 && num_LS == maxLSstep
        warning('gamma < 0')
    elseif num_LS == maxLSstep
            warning('gamma is too close to 0')        
	end    

    % continue the loop
    if (num_LS == maxLSstep) || updateS==1
         doLoop = 0;
    else
         num_LS = num_LS + 1; 
    end
end

if updateS~=1
	% use pure bfgs update if line search fails ?
	warning('cannot find telde alpha!')
    couldnotfindtildealpha=1;
end

if (num_LS == maxLSstep) && (gamma>-1e-6)
% 	sbfgsStep.size		= 0;   
% 	sbfgsIter   	  	= update_iterates(currIter,sbfgsStep,sbfgsProb,par);
% 	sk = 0;
% 	yk = 0;
% 	ak = 0;
% 	bk = 0;
% 	alpha = 0;
% 	gamma = 0;
     
    par.skipBFGS = 1;
	updateS=1;
% 	warning('skip current bfgs update \n');
end

wolfeProb.sk 	= sk;
wolfeProb.yk 	= yk;
wolfeProb.ak 	= ak;
wolfeProb.bk 	= bk;
wolfeProb.alpha = alpha;
wolfeProb.gamma = gamma;
    
ret = updateS;









