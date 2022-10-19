%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% test wolfe conditions
% 
function ret = test_wolfe(iter_new,step,prob_new,prob,par)

 c1 		= par.c1;
 c2 		= par.c2;

 obj_new	= prob_new.obj;
 grad_new	= prob_new.obj_grad;

 obj_old	= prob.obj;
 grad_old	= prob.obj_grad;

 step_size  = step.size;
 p			= step.x;


 % test 
 ret = 0;
 if obj_new <= obj_old + c1*step_size*grad_old'*p
	if grad_new'*p >= c2*grad_old'*p
	  ret=1;
	end
 end 

