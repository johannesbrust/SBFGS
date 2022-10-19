%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% set default option for S-BFGS
%
% 
function par = get_options()
% convergence
   	par.MaxIter 		= 1000;             	% max number of iterations (50)
    par.Tol 			= 1.e-6;			% max number of main loop tolerance (1.e-6)
    par.PrintLV         = 1;                % (1) to 10
    
% algorithm --- main
	par.alg       		= 'exact';    		% do exact, bfgs or s-bfgs: (bfgs), 
%     par.alg       		= 'bfgs';
    par.alg       		= 's-bfgs';    
   	par.addUnknown      = 'addNonLinObj';       % way to add unknown part: addQuadObj, addNonLinObj, addNonLinObj_Proj, (setRatio)
	par.exactRatio		= 0.25;				% ratio to add known part, valid if alg = 'setRatio'. (0.5)
    
	par.initB0			= 1;				% initial B0 = this_val * I. (1) 
	par.num_unknown     = 2;				% number of unknown variables, valid if alg = 'addQuadObj'/'addNonLinObj'. (2)
    
    par.initX0          = 0;                % initial point: (0):use ampl input, 1: set to one 

	par.checkInertia	= 1;				% check inertia or not: 0: not check (1)check 2: check descent direction

	par.init_reg 		= 1e-4;				% init primal regularizaion 
    
    par.safeguarding    = 0;                % apply safe guard in the zoom function

	par.sbfgsAlg       	= 2;				% way to build bfgs approximation: 1: A_k+K_k, (2): A_k+K_(k+1), only valid if alg='s-bfgs'

% about line search
	par.maxLSstep		= 30;

% about wolfe condition
	par.c1				= 0.0001;
	par.c2				= 0.9;

% do loop to increase exactRatio
    par.doLoop       	= 0;				% (0): not do; 1: do









% inner line-search steps
if par.maxLSstep == 1
par.innerFR = 1;
else
    for i = 0:par.maxLSstep-1
    par.innerFR(i+1) = 2^(-i);
    end
end

% inner line-search steps
par.modifyFR(1) = 0.9;
for i = 1:par.maxLSstep-1
    par.modifyFR(i+1) = par.modifyFR(i) + 0.9*10^(-i);
end


% last primal regularization
par.last_reg = 0;
% do inertia correction
par.done_correction = 0;

par.zoomLS = 0;

par.skipBFGS=0;

par.stopBFGS=0;

par.numFact = 0;

end



