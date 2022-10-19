%%%%%%% compute_step_CSBM1 (compact-structured-BFGS Minus Version 1) %%%%%%%%%
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. Structured objective functions are assumed to be
%
% f(x) = k(x) + u(x), x \in R^n,
%
% where k(x) has a known Hessian, i.e., k''(x) = K(x), and u(x) has an 'unknown'
% Hessian, i.e., u''(x) = U(x). The compact form of compact-structured-BFGS Minus
% (csBM) is
%
% B = Psi0 - Psi inv(M) Psi',
%
% where 
%
% Psi0  = A0 + K0 (typically A0 identity)
% Psi   = [Psi0*Sk, Yk]
% M     =   | Sk' Psi0 Sk,  Lk  |
%           | Lk',          -Dk |
% Sk'Yk = Lk + Rk
% Dk    = diag(Sk'Rk)
% Yk    = [y_m, ..., y_{k-1}]
% Sk    = [s_m, ..., s_{k-1}]
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leffeyer, M.Anitescu.
% A report on the compact representations is in DOCS/report_compact_sQN_011819 
%
% This function uses the interface of 'solvenlp.m' first written by
% Nai-Yuan. NOTE: This initial implementation is intended to compute steps using
% the compact formulas, alongside the 'original' recursive formulas.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [step,prob,par] = compute_step_CSBM1(nIter,prob_in,par_in,numIter,...
                            m,S,Q,R,D,YQ)

% Linear solve options
optsLSUTT.UT        = true; % Upper triangular
optsLSUTT.TRANSA    = true; % Transpose

optsLSUT.UT         = true; % Upper triangular

optsLSS.SYM         = true; % Symmetric

% get info
prob = prob_in;
par = par_in;
alg	= par.alg;
n   = prob.n;
	%m	= 0;

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
  
%   rhs = -prob.obj_grad;
% 
%   if par.checkInertia == 2
%     % check if descent     
%     [Mat_new,par] = check_inertia(n,0,Mat,rhs,prob,par,numIter);
%     step.x = Mat_new\rhs;
%   elseif par.checkInertia == 1 
%   % check and correct inertia        
%     [Mat_new,par] = check_inertia(n,0,Mat,rhs,prob,par,numIter);
%     step.x = Mat_new\rhs;
%   else
%     step.x = Mat\rhs;
%   end 

% Step computation
Psi0    = prob.Psi0;
g       = prob.obj_grad;

cidx    = nIter;

if nIter > m
    cidx = m;
end

buff                    = zeros(2*m,1);

p1                      = zeros(m,1);
p2                      = zeros(m,1);

p1(1:cidx)              = -(S(:,1:cidx)'*g);
p2(1:cidx)              = -(Q(:,1:cidx)'*g);
%buff(1:cidx)            = -(S(:,1:cidx)'*g);
%buff(cidx+1:2*cidx)     = -(Q(:,1:cidx)'*g);

buff(cidx+1:2*cidx)     = linsolve(R(1:cidx,1:cidx),p1(1:cidx),optsLSUT);
buff(1:cidx)            = linsolve(R(1:cidx,1:cidx),...
                        ((diag(D(1:cidx))+YQ(1:cidx,1:cidx))*buff(cidx+1:2*cidx)-p2(1:cidx)),optsLSUTT);


step.x                  = linsolve(Psi0,-g,optsLSS) + S(:,1:cidx)*buff(1:cidx) -...
                            Q(:,1:cidx)*buff(cidx+1:2*cidx);







