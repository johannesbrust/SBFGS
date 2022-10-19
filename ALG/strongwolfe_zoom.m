%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% line search method for wolfe condition
% 
%

function [zoomFlag,wolfe_step,par] = strongwolfe_zoom(iter,step,prob,par_in,func,alp_lo_in,alp_hi_in)

	alg			= par_in.alg;
    par         = par_in;

	trialIter	= iter;
    currIter	= iter;
	trialStep	= step;

 	obj_0		= prob.obj;	
	grad_0		= prob.obj_grad;
 	p			= trialStep.x;
	alp_lo		= alp_lo_in;
	alp_hi		= alp_hi_in;
   
	lo_Step		= step; 
	hi_Step		= step;  

    wolfe_step  = 0;
    
% alg parameters
	maxLSstep		= par.maxLSstep;
 	c1 				= par.c1;
 	c2 				= par.c2;
    safeG           = par.safeguarding;
	
% line-search info
    num_LS = 1;

% start line search loop  
	doLoop  = 1;
	maxflag = 0;
	ret=0;

	strong_wolfe_rhs = -c2*(grad_0'*p);
    
    % evaluate functions at the lo step and the hi step
	lo_Step.size	= alp_lo;
	lo_Iter   		= update_iterates(currIter,lo_Step,prob,par);
    lo_prob  		= evalFunc(lo_Iter,prob,par,func,1);
	obj_lo			= lo_prob.obj;
	gradp_lo		= lo_prob.obj_grad'*p;

	hi_Step.size	= alp_hi;
	hi_Iter   		= update_iterates(currIter,hi_Step,prob,par);
    hi_prob  		= evalFunc(hi_Iter,prob,par,func,1);
	obj_hi			= hi_prob.obj;   
    gradp_hi		= hi_prob.obj_grad'*p;
    
while doLoop == 1 

	% the minimizer of the quadratic interpolation
	aa = (obj_hi-obj_lo-(alp_hi-alp_lo)*gradp_lo)/((alp_hi-alp_lo)^2);
	bb = gradp_lo - 2*aa*alp_lo;

	alp_trial = -bb/(2*aa);
    
    % add safeguarding ???
    if safeG==1 && alp_trial < alp_lo + 0.1*abs(alp_hi-alp_lo)
        alp_trial = alp_lo + 0.1*abs(alp_hi-alp_lo);
        alp_lo;
        alp_hi;
    end
    
    % get trial iterate and evaluate functions
	trialStep.size	= alp_trial;    
	trialIter   	= update_iterates(currIter,trialStep,prob,par);
    trial_prob  	= evalFunc(trialIter,prob,par,func,1);

	% short-cut
 	obj_trial		= trial_prob.obj;
 	gradp_trial		= trial_prob.obj_grad'*p;

	% test strong Wolfe condition
 	if (obj_trial > obj_0 + c1*alp_trial*grad_0'*p) || (obj_trial > obj_lo)
		alp_hi   = alp_trial;
        obj_hi   = obj_trial;
        gradp_hi = gradp_trial;
	else
		if abs(gradp_trial) <= strong_wolfe_rhs
	  		wolfe_step = alp_trial;
			break;
		end
		if gradp_trial*(alp_hi - alp_lo) >= 0
			alp_hi   = alp_lo;
            obj_hi   = obj_lo;
            gradp_hi = gradp_lo;
		end 
		alp_lo   = alp_trial;
        obj_lo   = obj_trial;
        gradp_lo = gradp_trial;
 	end 
   
	if (num_LS == 10000)
       	maxflag = 1;
		break;
	end 
    num_LS 		= num_LS + 1;    
end

if num_LS > par.zoomLS
  par.zoomLS = num_LS;
end

zoomFlag=1;
if maxflag == 1
	% take full step if line search fails ?
	warning('strong wolfe zoom function fail. MAX ZOOM-IN ITER.');
	zoomFlag=0;
end
