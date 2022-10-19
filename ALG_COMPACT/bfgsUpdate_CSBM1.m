%%%%%%% bfgsUpdate_CSBM1 (compact-structured-BFGS Minus Version 1) %%%%%%%%%
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
% This function uses extends the interface of 'bfgsUpdate.m' first written by
% Nai-Yuan Chiang.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 04/03/19
function [ret,prob_new] = bfgsUpdate_CSBM1(newiter,iter,prob_new_in,prob,par,... % Original inputs
    m,nIter,S,Q,Y,Qtil,YQ,SQtil,R,L,D,Psi0)
    

%nIterp1         = nIter+1; % Shift for indexing

alg				= par.alg;
bfgsAlg			= par.sbfgsAlg;
Pmat            = prob.Pmat;

prob_new        = prob_new_in;

n               = prob.n;
% x_new           = newiter.x;
% x_old           = iter.x;
% grad_new		= prob_new.obj_grad_u;
% grad_old		= prob.obj_grad_u;
Mat			  	= prob.Bk;


ret             = 1;

gamma = 0;
if par.skipBFGS==0
    
    gamma  	= prob_new.gamma;
    
    if gamma <=0 && (bfgsAlg==1 && ~strcmp(par.addUnknown,'addNonLinObj_Proj') ) % for other cases, we can do regularization
        par.skipBFGS=1;
    end
    
end

% update B

if (par.stopBFGS==1)
    warning('set bfgs approximat matrix to 0');
    Mat_new = zeros(n);
else
    
    %   if (par.skipBFGS==1) && (  strcmp(alg,'bfgs') || ( strcmp(alg,'s-bfgs') && bfgsAlg==1)   )
    if (par.skipBFGS==1)
        
        warning('skip current bfgs update');
        
        ret=0;
        Mat_new = Mat;
        prob_new.Hes = prob.Hes;
        ret=1;
        
    else
        
        sk 		= prob_new.sk;
        ak 		= prob_new.ak;
        bk 		= prob_new.bk;
        alpha  	= prob_new.alpha;
        
        qk      = Psi0\bk;
        qktil   = Psi0*sk;
        
        cidx    = nIter+1;
        
        if cidx > m
            cidx = m;
        end
        
        %buff    = zeros(m,1);
        % Updates.
        if((nIter+1) <= m)
            
            S(:,cidx)           = sk;
            Q(:,cidx)           = qk;            
            Y(:,cidx)           = bk;
            Qtil(:,cidx)        = qktil;
            
        else
            
            S(:,1:cidx-1)           = S(:,2:cidx);
            Q(:,1:cidx-1)           = Q(:,2:cidx);
            R(:,1:cidx-1)           = R(:,2:cidx);
            D(1:cidx-1)             = D(2:cidx);
            YQ(1:cidx-1,1:cidx-1)   = YQ(2:cidx,2:cidx);
            
            L(1:cidx-1,1:cidx-1)    = L(2:cidx,2:cidx);
            SQtil(1:cidx-1,1:cidx-1)= SQtil(2:cidx,2:cidx);
            
            S(:,cidx)               = sk;
            Q(:,cidx)               = qk;
            Y(:,cidx)               = bk;
            Qtil(:,cidx)            = qktil;
            
        end
        
        R(1:cidx,cidx)      = S(:,1:cidx)'*bk; %buff(1:cidx);
        D(cidx,1)           = R(cidx,cidx);
        
        YQ(cidx,1:cidx)     = bk'*Q(:,1:cidx);
        YQ(1:cidx,cidx)     = YQ(cidx,1:cidx)';
        
        L(cidx,1:cidx-1)    = sk'*Y(:,1:cidx-1);
        
        SQtil(cidx,1:cidx)  = sk'*Qtil(:,1:cidx);
        SQtil(1:cidx,cidx)  = SQtil(cidx,1:cidx)';
        
        %Mat_new = Mat - (ak*ak')/alpha + bk*bk'/gamma;
        
        if bfgsAlg==1 && strcmp(alg,'s-bfgs')
            if strcmp(par.addUnknown,'addNonLinObj_Proj')
                Mat_new = Mat_new + Pmat*prob.Hes*Pmat' - Pmat*prob_new.Hes*Pmat';
            else
                
                Mat_new = Psi0 - prob_new.Hes - ...
                            [Qtil(:,1:cidx),Y(:,1:cidx)] * ...
                            (([SQtil(1:cidx,1:cidx),L(1:cidx,1:cidx);...
                                L(1:cidx,1:cidx)',-diag(D(1:cidx,1))])\[Qtil(:,1:cidx),Y(:,1:cidx)]');
                            
                %Mat_new = Mat_new + prob.Hes - prob_new.Hes;
                
            end
        end
    end
    
end

prob_new.Bk = Mat_new;
prob_new.S      = S;
prob_new.Q      = Q;
prob_new.R      = R;
prob_new.D      = D;
prob_new.YQ     = YQ;
prob_new.Y      = Y;
prob_new.Qtil   = Qtil;
prob_new.L      = L;
prob_new.SQtil  = SQtil;



% check if new hessian approximation is pd
%  if strcmp(alg,'bfgs')
%    Hes_use = zeros(dim_x,dim_x);
%  elseif bfgsAlg==1 && strcmp(alg,'s-bfgs')
%   	Hes_use   = prob.Hes;
%  elseif bfgsAlg==2 && strcmp(alg,'s-bfgs')
%	Hes_use   = prob_new.Hes;
%  else
%   	error('bfgsAlg should be 1 or 2 within s-bfgs ');
%  end
%
%   oldAppro  = Mat+Hes_use;
%   TestMat1  = (eye(dim_x)-sk*bk'/gamma);
%   TestMat   = TestMat1*(inv(oldAppro))*TestMat1'+(sk*sk')/gamma;
%   newAppro1 = inv(TestMat);
%
%   if ~strcmp(alg,'bfgs')
%     newAppro2 = Mat_new + prob_new.Hes;
%   else
%     newAppro2 = Mat_new;
%   end
%
%   checek = 11;
%   if  ~all(eig(oldAppro)>0)
%       warning('Current Hessian approximation is not PD');
%   end
%   if  ~all(eig(newAppro2)>0)
%       warning('New Hessian approximation is not PD');
%   end
%
%   if max(max(abs(newAppro1-newAppro2))) >1e-3
%      warning('nemerical error! cond(mat1) = %2.2e, cond(mat2) = %2.2e, difference: %3.6e', cond(newAppro1),cond(newAppro2),max(max(abs(newAppro1-newAppro2))) );
%   end






