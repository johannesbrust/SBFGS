%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test wolfe conditions
% 
function interXMin = minInterpolation(alp_lo,alp_hi,lo_prob,hi_prob,iter_in,step_in,prob_in,par,func)

	trialStep	= step_in;

	p			= trialStep.x;

	obj_lo			= lo_prob.obj;
	obj_hi			= hi_prob.obj;

	gradTimesP_lo	= (lo_prob.obj_grad)' * p;
	gradTimesP_hi	= (hi_prob.obj_grad)' * p;

	d1				= gradTimesP_lo + gradTimesP_hi - 3 *(obj_lo-obj_hi)/(alp_lo-alp_hi);
	signD			= 1;
	if alp_hi-alp_lo < 0
		signD = -1;
	end
	d2				= signD*sqrt( d1^2 - gradTimesP_lo * gradTimesP_hi );

	alp_min			= alp_hi-(alp_hi-alp_lo)*( (gradTimesP_hi+d2-d1)/(gradTimesP_hi-gradTimesP_lo+2*d2)   );


	interXMin 		= 	alp_lo;
	if obj_hi < obj_lo
        interXMin 		= 	alp_hi;             
    end
    
    if  alp_min>0
        trialStep.size	= alp_min;
        trialIter   	= update_iterates(iter_in,trialStep,prob_in,par);
        trial_prob  	= evalFunc(trialIter,prob_in,par,func,1);
        obj_trial		= trial_prob.obj;
        
        if obj_trial < obj_hi && obj_trial < obj_lo
            interXMin 		= 	alp_min;
        end

    end
