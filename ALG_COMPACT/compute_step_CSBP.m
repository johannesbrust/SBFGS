%%%%%%% compute_step_CSBP (compact-structured-BFGS Plus Version 1 Dense) %%
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. Structured objective functions are assumed to be
%
% f(x) = k(x) + u(x), x \in R^n,
%
% where k(x) has a known Hessian, i.e., k''(x) = K(x), and u(x) has an 'unknown'
% Hessian, i.e., u''(x) = U(x). The compact form of compact-structured-BFGS
% Plus (CSBM) is
%
% B = (K_{k+1} + A0) - Psi inv(M) Psi',
%
% where 
%
% Psi   = [Qk, Yk]
% Qk    = Vk + A0*Sk
% M     =   | Dkv + Lkv + Lkv' + Sk' A0 Sk,     Lk  |
%           | Lk',                              -Dk |
% Sk'Yk = Lk + Rk
% Dk    = diag(Sk'Yk)
% Sk'Vk = Lkv + Rkv
% Dkv   = diag(Sk'Vk) 
% Yk    = [y_m, ..., y_{k-1}]
% Sk    = [s_m, ..., s_{k-1}]
% Vk    = [K_{m+1}s_m, ..., K_ks_{k-1}]
%
%
% Modification of the 'compute_step' function.
%
% Initial contributors: J.J.Brust, C.G.Petra, S. Leyffer, M.Anitescu.
% A report on the compact representations is in DOCS/report_compact_sQN_011819 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [step,prob,par] = compute_step_CSBP(nIter,prob_in,par_in,numIter)

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
    [Mat_new,par] = check_inertia_CSBP(n,0,Mat,rhs,prob,par,numIter);
    step.x = Mat_new\rhs;
  else
    step.x = Mat\rhs;
  end 

  



