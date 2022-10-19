%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The optimization problem:
% min   f(x) = fk(x) + fu(x)
%
% 
function [bool_conv,sol] = solvenlp(prob,par,func)

%% Extra outputs (J.B., 04/18/19)
ctime               = tic;  % Timer

dim_x               = prob.n;

iter.x              = prob.x0;

%%% get algorithmic options
maxIter		= par.MaxIter;
OutTol		= par.Tol;
alg			= par.alg;

%%% eval functions & check init point
prob = evalFunc(iter,prob,par,func,2);
if ~strcmp(alg,'exact') 
  prob = initB0(prob,par);
end

doMainLoop = 1;
residual = norm(prob.obj_grad,Inf);
if residual <= OutTol
  doMainLoop = 0;
  num_LS     = 0; 
end

%%% start main loop
ret = 1;
nIter 		= 0;
if(par.PrintLV>1)
    fprintf('Iter       obj          Resid        ||d||         alpha       #LS UseZoom Reg     PriR\n');
    fprintf('%d\t%e\t  %.2e	\n', nIter, prob.obj, residual);    
end
while (nIter < maxIter) && (doMainLoop == 1)

    %%% solve linear system
    [steps,prob,par] = compute_step(nIter,prob,par,nIter);
    
	%%% do line search
%     [ret, newIter,newStep,prob_new,num_LS] = line_search(iter,steps,prob,par,func);
%     [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_strongwolfe(iter,steps,prob,par,func);
    [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_jorge(iter,steps,prob,par,func);
    
    if ret~=1
        ctime = toc(ctime);
       break; 
    end
    
	if ~strcmp(alg,'exact') 
  	  [ret,prob_new] = bfgsUpdate(newIter,iter,prob_new,prob,par);
      if ret~=1
          ctime = toc(ctime);
        break; 
      end
    end
    
	%%% check residuals
	residual = norm(prob_new.obj_grad,Inf);
	if residual <= OutTol
      doMainLoop = 0;
    end

	nIter = nIter+1; 
    
    if(par.PrintLV>1)
        if par.done_correction ~=0
            priReg=log10(par.last_reg);
            fprintf('%d\t%e\t  %.2e\t%.2e\t%.2e\t%d\t%d\t%d\t%.1f\n', nIter, prob_new.obj, residual, norm(newStep.x,Inf),newStep.size,num_LS,use_zoom,par.done_correction,priReg);        
        else
            priReg='-';
            fprintf('%d\t%e\t  %.2e\t%.2e\t%.2e\t%d\t%d\t%d\t%s\n', nIter, prob_new.obj, residual, norm(newStep.x,Inf),newStep.size,num_LS,use_zoom,par.done_correction,priReg);
        end
    end

    prob  = prob_new;
    iter  = newIter;
    steps = newStep;

end

if doMainLoop == 0
	bool_conv = 1;
else
	bool_conv = 0;
end

sol = iter;
sol.nIter = nIter;
sol.name = 'probname';
sol.zoomLS = par.zoomLS;
sol.obj = prob.obj;
sol.numFact = par.numFact;
% sol.num_LS = num_LS;
if (bool_conv == 1); sol.ctime = toc(ctime); else sol.ctime = ctime; end

sol.OK = 0;
if bool_conv==1
   sol.flag = 'OK  ';
   sol.OK = 1;
elseif nIter == maxIter
   sol.flag = 'MAX '; 
   sol.OK = 2;
else
   sol.flag = 'Fail'; 
   sol.OK = 0;
end

