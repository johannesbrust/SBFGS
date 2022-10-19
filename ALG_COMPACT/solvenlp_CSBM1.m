%%%%%%% solvenlp_CSBM1 (compact-structured-BFGS Minus Version 1) %%%%%%%%%
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
% Initial version: J.B., 04/03/19

function [bool_conv,sol] = solvenlp_CSBM1(prob,par,func)

%dim_x       	   = prob.n;
iter.x             = prob.x0;

%%% get algorithmic options
maxIter		= par.MaxIter;
OutTol		= par.Tol;
alg			= par.alg;

%%% eval functions & check init point
prob = evalFunc(iter,prob,par,func,2);
if ~strcmp(alg,'exact') 
  prob = initB0(prob,par);
  
  % Store initial matrix Psi0
  Psi0 = prob.Bk + prob.Hes;
  
end

doMainLoop  = 1;
g           = prob.obj_grad;
residual    = norm(g,Inf);
if residual <= OutTol
  doMainLoop = 0;
  num_LS     = 0; 
end

% Initializations
% Preparations for inverse computations
m           = par.m; % Memory parameter
n           = prob.n;
S           = zeros(n,m);
Q           = zeros(n,m); % inv(Psi0)Y
R           = zeros(m,m);
D           = zeros(m,1);
YQ          = zeros(m,m); % Y'*Q = Y'*inv(Psi0)*Y 

% Matrices for original representations
Y           = zeros(n,m);
Qtil        = zeros(n,m); % Psi0*S
L           = zeros(m,m);
SQtil       = zeros(m,m); % S'*Qtil = S'*Psi0*S

% Errors
eBs         = zeros(m,1); % Errors recursive and compact forms (matrices)
eSs         = zeros(m,1); % Errors steps (vectors)
%Psi0g       = Psi0\g;


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
    
    %%% Solve linear system compact representation
    [stepsc,probc,parc] = compute_step_CSBM1(nIter,prob,par,nIter,...
                            m,S,Q,R,D,YQ);
    
    % Step equivalence check
    if((nIter+1) <= m)               
        
        eS              = norm(steps.x-stepsc.x);
        eSs(nIter+1)    = eS;
        
    end
    
                        
	%%% do line search
%     [ret, newIter,newStep,prob_new,num_LS] = line_search(iter,steps,prob,par,func);
%     [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_strongwolfe(iter,steps,prob,par,func);
    [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_jorge(iter,steps,prob,par,func);
    
    if ret~=1
       break; 
    end
    
	if ~strcmp(alg,'exact') 
        
       % Recursive update. Also updates Psi0.
  	  [ret,prob_new] = bfgsUpdate(newIter,iter,prob_new,prob,par);
      
      % Compact representation
      [retc,prob_newc] = bfgsUpdate_CSBM1(newIter,iter,prob_new,prob,par,... 
                m,nIter,S,Q,Y,Qtil,YQ,SQtil,R,L,D,prob.Psi0);
            
       S    = prob_newc.S;
       Q    = prob_newc.Q;
       YQ   = prob_newc.YQ;
       R    = prob_newc.R;
       D    = prob_newc.D;
       Y    = prob_newc.Y;
       Qtil = prob_newc.Qtil;
       L    = prob_newc.L;
       SQtil= prob_newc.SQtil;
       
       prob_new.Psi0 = prob.Psi0;
              
       % Matrix equivalence check
       if ((nIter+1) <= m)
           
            eB              = norm(prob_new.Bk-prob_newc.Bk,'fro');
            eBs(nIter+1)    = eB;
            
       end
       
      if ret~=1
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

