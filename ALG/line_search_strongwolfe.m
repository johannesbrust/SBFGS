%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% line search method for wolfe condition
% 
%

function [ret,newIter,newStep,newProb,num_LS,use_zoom,par] = line_search_strongwolfe(iter,step,prob,par_in,func)

	alg			= par_in.alg;
    par         = par_in;
	bfgsAlg		= par.sbfgsAlg;

	trialIter	= iter;
    currIter	= iter;
	trialStep	= step;

	alp_max   	= 2.0;
	alp_init   	= 1.0;

 	obj_0		= prob.obj;	
	grad_0		= prob.obj_grad;
 	p			= trialStep.x;
    
% alg parameters
	maxLSstep		= par.maxLSstep;
 	c1 				= par.c1;
 	c2 				= par.c2;

% alpha_last is alpha_{i-1}
% alpha_curr is alpha_i
	alp_last 	= 0;
	alp_trial  	= alp_init;

% line-search info
    num_LS = 1;

% check if direction is descent
    if grad_0'*p >=0
        aws=grad_0'*p
        warning('Not decsent direction');
    end
    
% start line search loop  
	doLoop  = 1;
	maxflag = 0;
	findwolfe=0;
    zoomFlag = 1;
    ret2 = 1;

	strong_wolfe_rhs = -c2*(grad_0'*p);
	obj_last = obj_0;
    
    use_zoom = 0;

while doLoop == 1 
    % get trial iterate and evaluate functions
	trialStep.size	= alp_trial;
	trialIter   	= update_iterates(currIter,trialStep,prob,par);
    trial_prob  	= evalFunc(trialIter,prob,par,func,1);

	% short-cut
 	obj_trial		= trial_prob.obj;
 	grad_trial		= trial_prob.obj_grad;

	% test strong Wolfe condition
	findwolfe = 0;
    if (obj_trial > obj_0 + c1*alp_trial*grad_0'*p) || (num_LS>1 && obj_trial >= obj_last)
    	[zoomFlag,trialStep.size,par] = strongwolfe_zoom(currIter,trialStep,prob,par,func,alp_last,alp_trial);
    	findwolfe = 1;
        use_zoom  = 1;
    end
	if findwolfe==0 && abs(grad_trial'*p) <= strong_wolfe_rhs
		findwolfe = 1;
        use_zoom  = 0;
	end  
	if findwolfe==0 && (grad_trial'*p) >= 0
		[zoomFlag,trialStep.size,par] = strongwolfe_zoom(currIter,trialStep,prob,par,func,alp_trial,alp_last);
	    findwolfe = 1;
        use_zoom  = 1;
	end

    % continue the loop
	if (findwolfe==1) 
	  doLoop = 0;
	  trialIter   		= update_iterates(currIter,trialStep,prob,par);
      trial_prob  		= evalFunc(trialIter,prob,par,func,2);
      if strcmp(alg,'s-bfgs')  && zoomFlag==1
        % find correct step-size for s-bfgs method
	    [ret2,trial_prob] = modify_WolfeStep(iter,prob,trialIter,trialStep,trial_prob,par,func);
      end
	  break;
    end

	if (num_LS == maxLSstep)
       doLoop = 0;
	   maxflag = 1;
	   break;
    end

	alp_last 	= alp_trial;
	obj_last	= obj_trial;
	alp_trial 	= (alp_max + alp_last)/2; 
    num_LS 		= num_LS + 1; 
end

newIter    	= trialIter;
newProb		= trial_prob;
newStep     = trialStep;

if maxflag == 1
	% take full step if line search fails ?
	warning('Wolfe condition fails. MAX LS ITER')
end

if findwolfe==1 && ret2==1 && zoomFlag==1
    ret = 1;
else
    ret=0;
end


