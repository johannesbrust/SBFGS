% compute step using eigenvalue inertia detection
%
% Nai-Yuan Chiang, 2015
%
function [step,prob,par] = compute_step(nIter,prob_in,par_in,numIter)

% get info
    prob = prob_in;
    par = par_in;
	alg	= par.alg;
    n   = prob.n;
	m	= 0;

% set iterate matrix
if strcmp(alg,'exact')
  Mat = prob.Hes;
elseif strcmp(alg,'bfgs')
  Mat = prob.Bk;
elseif strcmp(alg,'s-bfgs')
    
  if strcmp(par.addUnknown,'addNonLinObj_Proj')
      Mat = prob.Hes + prob.Pmat'*prob.Bk*prob.Pmat; 
  else
      Mat = prob.Hes + prob.Bk;
  end
 
  if(nIter==0)
    scalB0  = 10;
    trialNum = 0;      
    while ~all(eig(Mat)>0)
      trialNum = trialNum+1;
      Mat = prob.Hes + power(scalB0,trialNum)*eye(n);
%       all(scalB0>0);
%         scalB0 = abs(diag(prob.Hes)) + 1e-6*ones(n,1);
%         Mat =diag(scalB0);
%         prob.Bk = Mat;
%         prob.Hes= zeros(n);          
    end
    
    prob.Bk = power(scalB0,trialNum)*prob.Bk;
    
    % Update initial matrix
    prob.Psi0 = Mat;
    
  end
end
  
  rhs = -prob.obj_grad;

  if par.checkInertia == 2
    % check if descent     
    [Mat_new,par] = check_inertia(n,0,Mat,rhs,prob,par,numIter);
    step.x = Mat_new\rhs;
  elseif par.checkInertia == 1 
  % check and correct inertia        
    [Mat_new,par] = check_inertia(n,0,Mat,rhs,prob,par,numIter);
    step.x = Mat_new\rhs;
  else
    step.x = Mat\rhs;
  end 

  



