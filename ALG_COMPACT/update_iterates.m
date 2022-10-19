%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% update iterates
%
% 
function [iter_new] = update_iterates(iter,step,prob,par)

	iter_new     	= iter;
	iter_new.x      = iter.x   + step.size*step.x;









