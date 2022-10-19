%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% line search method 
% 
%

function [ret,newIter,newStep,newProb,num_LS] = line_search(iter,step,prob,par,func)

	alg			= par.alg;

	trialIter	= iter;
    currIter	= iter;
	trialStep	= step;

	alp_max   	= 1.0;
    
% alg parameters
	FR					= par.innerFR;
	maxLSstep			= par.maxLSstep;

% line-search info
    num_LS = 1;

% start line search loop  
	doLoop = 1;
	ret=0;
    ret2=1;
while doLoop == 1 
    % get trial iterate and merit 
    trialStep.size    = alp_max*FR(num_LS);   
	trialIter   	  = update_iterates(currIter,trialStep,prob,par);

    % compute derivatives
    trial_prob  = evalFunc(trialIter,prob,par,func,1);

	% test Wolfe condition
	ret = test_wolfe(trialIter,trialStep,trial_prob,prob,par);
        

    % continue the loop
	if (ret==1) 
	  doLoop = 0;
      if strcmp(alg,'s-bfgs') 
        % find correct step-size for s-bfgs method
	    [ret2,trial_prob] = modify_WolfeStep(iter,prob,trialIter,trialStep,trial_prob,par,func);
        if ret2 ~= 1
           break; 
        end
      end
    elseif (num_LS == maxLSstep)
       doLoop = 0;
	else
      num_LS = num_LS + 1; 
    end
end


	newIter    	= trialIter;
	newProb		= trial_prob;
    newStep     = trialStep;

if ret ~= 1
	% take full step if line search fails ?
	warning('Wolfe condition fails.')
end

if ret==1 && ret2==1
    ret = 1;
else
    ret=0;
end





