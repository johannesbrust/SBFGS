%%%%%%% CSBP2D_SCALE (compact-structured-BFGS Plus Version 2 Dense %%%%%%%%
%                       with Scaling)
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
% Version 2 (V2) sets 
%
% A0 = gamma_k.I = Non-constant. (I.e., A0 is a non-constant multiple of identity
% matrix)
%
% The scaling choices are:
%
% (1): delta = bk'*bk/sk'*bk,
% (2): delta = sk'*bk/sk'*sk,
% (3): delta = yk'*yk/sk'*yk,
% (4): delta = sk'*yk/sk'*sk,
%
% where bk = yk + K(xk1)sk, yk = grad(u(xk1)) - grad(u(xk)),
% and K(x) is the known Hessian, and the Hessian of u(x) is unknown.
%
% Initial contributors: J.J.Brust, C.G.Petra, S. Leyffer, M.Anitescu.
% A report on the compact representations is in DOCS/report_compact_sQN_011819 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initial version: J.B., 05/07/19
% 10/22/19, J.B., Version with scaling choices
% 11/14/19, J.B., Extension to store scaling data

function [bool_conv,sol] = CSBP2_D_SCALE(prob,par,func)

%% Extra outputs
ctime       = tic;  % Timer
%nhes        = 0;    % No. Hessian evaluations

iter.x      = prob.x0;

%%% get algorithmic options
maxIter		= par.MaxIter;
OutTol		= par.Tol;

alg         = par.alg;

%%% eval functions & check init point
prob        = evalFunc(iter,prob,par,func,2);
%nhes        = nhes + 1;



% 04/09/19, J.B., Sparse matrices
%spB         = sparse(prob.Bk);
%spK0        = sparse(prob.Hes);

%Psi0        = spB + spK0;
%CSB         = Psi0;

doMainLoop  = 1;
g           = prob.obj_grad;

residual    = norm(g,Inf);
if residual <= OutTol
  doMainLoop = 0;
end

%% Initializations
% Matrices to compute compact representations
m           = par.m; % Memory parameter
n           = prob.n;
S           = zeros(n,m);
Q           = zeros(n,m); % Vk + A0*Sk
%R           = zeros(m,m);
D           = zeros(m,1);
%YQ          = zeros(m,m); % Y'*Q = Y'*inv(Psi0)*Y 

Y           = zeros(n,m);
V           = zeros(n,m); % [ ..., K_{i+1}*si, ... ]
L           = zeros(m,m);
%SQtil       = zeros(m,m); % S'*Qtil = S'*Psi0*S
SS          = zeros(m,m); % S'*S
Dv          = zeros(m,1);
Lv          = zeros(m,m);

% Temporary storage buffers for step computations
buff        = zeros(2*m,1);
p1          = zeros(m,1);
p2          = zeros(m,1);

% Matrices to compute with B

midx        = zeros(m,1); % Memory indices, e.g. S(:,midx)
midx(1:m)   = 1:m;

%% Linear Solve Options
optsLSUTT.UT        = true; % Upper triangular
optsLSUTT.TRANSA    = true; % Transpose

optsLSUT.UT         = true; % Upper triangular

optsLSS.SYM         = true; % Symmetric

%% Line-search Options
optsLINESRCH.c1         = par.c1;
optsLINESRCH.c2         = par.c2;

optsLINESRCH.alp_max    = 2.0;
optsLINESRCH.alp_init   = 1.0;
optsLINESRCH.xtol       = 1e-6;
optsLINESRCH.stpmin     = 0;
optsLINESRCH.stpmax     = optsLINESRCH.alp_max;
optsLINESRCH.maxfev     = 3000;

%% Storing Options
deltas          = [];
storeDeltas     = 0;
if isfield(par,'storeDeltas')
    if par.storeDeltas == true
        storeDeltas   = par.storeDeltas;
        deltas        = zeros((maxIter+1),1);
    else
        storeDeltas   = false;
    end
end

%% Eigenvalue check and initialization
% NonLinObj_Proj currently not supported
% if strcmp(par.addUnknown,'addNonLinObj_Proj')
%     
%     CSB = prob.Hes + prob.Pmat'*prob.Bk*prob.Pmat;
%     
% end

% 04/26/19, J.B., Dense matrices
% Multiple of identity initial matrix
%prob        = initB0(prob,par); 

In          = eye(n);
On          = ones(n,1);
gaminit     = par.initB0;
B0          = gaminit * In;

K0          = prob.Hes;
Psi0        = B0 + K0;

prob.last_obj_grad_u = 0;
prob.last_obj_grad_k = 0;

%spIn        = speye(n);

scalB0      = 10;
trialNum    = 0;

lam1        = eigs(Psi0,1,'SA'); % Smallest algebraic eigenvalues

while lam1 < 0 

    %while ~all(eig(CSB)>0) 
    
    
    trialNum    = trialNum+1;
    gaminit     = power(scalB0,trialNum);
    B0          = gaminit*In;
    Psi0        = K0 + B0;
    
    lam1        = eigs(Psi0,1,'SA');
    
    %       all(scalB0>0);
    %         scalB0 = abs(diag(prob.Hes)) + 1e-6*ones(n,1);
    %         Mat =diag(scalB0);
    %         prob.Bk = Mat;
    %         prob.Hes= zeros(n);
end

prob.Bk     = B0;

%prob.Bk     = power(scalB0,trialNum)*spB;
%prob.Bk     = power(scalB0,trialNum)*full(spB)+full(spK0);
%prob.Bk     = full(CSB);

%prob.Bk     = full(spIn);

%prob.Psi0   = CSB;
%Psi0        = CSB;
%fPsi0       = full(Psi0);
%Psi0        = full(CSB);

%%% start main loop
ret         = 1;
nIter 		= 0;

if(par.PrintLV>1)
    fprintf('Iter       obj          Resid        ||d||         alpha       #LS UseZoom Reg     PriR\n');
    fprintf('%d\t%e\t  %.2e	\n', nIter, prob.obj, residual);    
end

% 05/06/19, J.B.
% Initially compute original steps next to compact steps.

while (nIter < maxIter) && (doMainLoop == 1)

    %%% solve linear system
    [steps,prob,par] = compute_step(nIter,prob,par,nIter);
    %[steps,prob,par] = compute_step_CSBP(nIter,prob,par,nIter);
    
    %% Step computation
%     cidx    = nIter;
% 
%     if nIter > m
%         cidx = m;
%     end

    % REQUIRE
    % Limited-memory matrices
    % g
    
%     p1(1:cidx)              = -(S(:,midx(1:cidx))'*g);
%     p2(1:cidx)              = -(Q(:,midx(1:cidx))'*g);
% 
%     %buff(cidx+1:2*cidx)     = linsolve(R(1:cidx,1:cidx),p1(1:cidx),optsLSUT);
%     
%     buff(cidx+1:2*cidx)     = R(1:cidx,1:cidx)\p1(1:cidx);
%     
% %     buff(1:cidx)            = linsolve(R(1:cidx,1:cidx),...
% %                             ((diag(D(1:cidx))+YQ(1:cidx,1:cidx))*buff(cidx+1:2*cidx)-p2(1:cidx)),optsLSUTT);
% 
%     buff(1:cidx)            = R(1:cidx,1:cidx)'\((diag(D(1:cidx))+YQ(1:cidx,1:cidx))*buff(cidx+1:2*cidx)-p2(1:cidx));
%                             
% 
% 
% %     st                      = Psi0 \ (-g) + ... % Sparse 
% %                                 S(:,midx(1:cidx))*buff(1:cidx) -...
% %                                 Q(:,midx(1:cidx))*buff(cidx+1:2*cidx);
% 
% %     st                      = linsolve(Psi0,-g,optsLSS) + ...  
% %                                 S(:,midx(1:cidx))*buff(1:cidx) -...
% %                                 Q(:,midx(1:cidx))*buff(cidx+1:2*cidx);
% 
% %     st                      = full(Psi0) \ (-g) + ... % full 
% %                                 S(:,midx(1:cidx))*buff(1:cidx) -...
% %                                 Q(:,midx(1:cidx))*buff(cidx+1:2*cidx);
% 
%     st                      = Psi0 \ (-g) + ... % full 
%                                 S(:,midx(1:cidx))*buff(1:cidx) -...
%                                 Q(:,midx(1:cidx))*buff(cidx+1:2*cidx);
% 
%     steps.x                 = st;                       
    
    %steps.x                  = stepst.x;
    
    % Error 
    %e = norm(stepst.x-st)/min(norm(st),norm(stepst.x));
    
	%%% do line search
%     [ret, newIter,newStep,prob_new,num_LS] = line_search(iter,steps,prob,par,func);
%     [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_strongwolfe(iter,steps,prob,par,func);
    
    % REQUIRE
    % Line-search parameters
    % st, Step
    % Iterate

    [ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_jorge(iter,steps,prob,par,func);
    %[ret, newIter,newStep,prob_new,num_LS,use_zoom,par] = line_search_CSBM1(iter,steps,prob,par,func,optsLINESRCH);
    
    %nhes = nhes + prob_new.nhes;
    
    if ret~=1
        
        ctime = toc(ctime);
        
       break; 
    end
    
    % REQUIRE
    % bk
    % sk
    % Hes
    
    %% Limited-Memory BFGS update
    
    bfgsAlg			= par.sbfgsAlg;
    Pmat            = prob.Pmat;
       
    Mat			  	= prob.Bk;
    
    ret             = 1;
    
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
        
        if (par.skipBFGS==1) && (  strcmp(alg,'bfgs') || ( strcmp(alg,'s-bfgs') && bfgsAlg==1)   )
        %if (par.skipBFGS==1)
            
            warning('skip current bfgs update');
            
            Mat_new = Mat;
            prob_new.Hes = prob.Hes;
            ret=1;
            
        else
            
            sk 		= prob_new.sk;
            bk 		= prob_new.bk;
            
            vk      = prob_new.Hes*sk;
            
            % 05/07/19, J.B., Gamma scaling 
%             gammfac  = abs(gamma);
%             %gammfac = gamma;
%             epsgam  = 1.0e+11;
%             
%             gammfac = min(gammfac,epsgam);
%             gammfac = max(gammfac,1/epsgam);
            
            % 10/22/19, J.B. switch over scaling choices
            
            gaminit = (bk'*bk)/gamma;
            
            switch par.Scaling
                case 2
                    yk      = prob_new.yk;
                    gamma   = sk'*yk;
                    gaminit = (yk'*yk)/gamma;                    
                case 3
                    ss      = sk'*sk;
                    gaminit = gamma/ss;
                case 4
                    yk      = prob_new.yk;
                    gamma   = sk'*yk;
                    gaminit = gamma/(sk'*sk);
            end
            
            if storeDeltas == true
                deltas(nIter+2,1) = gaminit;
            end
            
%             if agamma > epsgam
%                
%                 gammfac = epsgam;
%                 
%             elseif agamma < 1/epsgam
%                 
%                 gammfac = 1/epsgam;
%                 
%             end
            
            
            
            %qk      = Psi0\bk;
            %qktil   = Psi0*sk;
            
            cidx    = nIter+1;
            
            if cidx > m
                cidx = m;
            end
            
            %buff    = zeros(m,1);
            % Updates.
            if((nIter+1) <= m)
                
                S(:,midx(cidx))           = sk;
                
                Y(:,midx(cidx))           = bk;
                V(:,midx(cidx))           = vk;
                
                
                %Qtil(:,midx(cidx))        = qktil;
                
                % For memory index test
                %S_ = S;
                
            else
                
                % Large arrays use 'midx' to avoid excess copying
                midx1                   = midx(1);
                midx(1:m-1)             = midx(2:m);
                midx(m)                 = midx1;
                
                S(:,midx1)              = sk;
                %Q(:,midx1)              = vk + gaminit.*sk;
                Y(:,midx1)              = bk;
                V(:,midx1)              = vk;
                %Qtil(:,midx1)           = qktil;
                
                %S(:,1:cidx-1)           = S(:,2:cidx);
                %Q(:,1:cidx-1)           = Q(:,2:cidx);
                
                % Small arrays 
                %R(1:cidx-1,1:cidx-1)    = R(2:cidx,2:cidx);
                D(1:cidx-1)             = D(2:cidx);
                
                L(1:cidx-1,1:cidx-1)    = L(2:cidx,2:cidx);                
                SS(1:cidx-1,1:cidx-1)   = SS(2:cidx,2:cidx);
                
                % V arrays
                Lv(1:cidx-1,1:cidx-1)   = Lv(2:cidx,2:cidx);
                Dv(1:cidx-1,1)          = Dv(2:cidx,1);
                
                % Memory index test
                
%                 S_(:,1:cidx-1)              = S_(:,2:cidx);
%                 S_(:,cidx)                  = sk;
%                 
%                 eS = norm(S(:,midx(1:cidx))-S_(:,1:cidx),'fro');
                
            end
            
            % Potentially large memory copy
            Q(:,1:cidx)         = V(:,midx(1:cidx)) + gaminit.*S(:,midx(1:cidx));
            
            %R(1:cidx,cidx)      = S(:,midx(1:cidx))'*bk; %buff(1:cidx);
            sb                  = sk'*bk;
            sv                  = sk'*vk;
            
            D(cidx,1)           = sb;
            Dv(cidx,1)          = sv;
            
            SS(cidx,1:cidx)     = (sk'*S(:,midx(1:cidx)));
            SS(1:cidx,cidx)     = SS(cidx,1:cidx)';
            
            L(cidx,1:cidx-1)    = sk'*Y(:,midx(1:cidx-1));
            
            Lv(cidx,1:cidx-1)   = sk'*V(:,midx(1:cidx-1));
            
            %SQtil(cidx,1:cidx)  = sk'*Qtil(:,midx(1:cidx));
            %SQtil(1:cidx,cidx)  = SQtil(cidx,1:cidx)';
            
            %Mat_new = Mat - (ak*ak')/alpha + bk*bk'/gamma;
            
            if strcmp(par.addUnknown,'addNonLinObj_Proj')
                Mat_new = Mat_new + Pmat*prob.Hes*Pmat' - Pmat*prob_new.Hes*Pmat';
            else
                
%                 Mat_new = Psi0 - prob_new.Hes - ...
%                     [Qtil(:,midx(1:cidx)),Y(:,midx(1:cidx))] * ...
%                     (([SQtil(1:cidx,1:cidx),L(1:cidx,1:cidx);...
%                     L(1:cidx,1:cidx)',-diag(D(1:cidx,1))])\[Qtil(:,midx(1:cidx)),Y(:,midx(1:cidx))]');
                % prob_new.Hes +
                Mat_new =  diag(gaminit.*On) - [Q(:,1:cidx),Y(:,midx(1:cidx))] * ...
                     (([Lv(1:cidx,1:cidx)+Lv(1:cidx,1:cidx)' + diag(Dv(1:cidx)) + gaminit.*SS(1:cidx,1:cidx),L(1:cidx,1:cidx);...
                     L(1:cidx,1:cidx)',-diag(D(1:cidx,1))])\[Q(:,1:cidx),Y(:,midx(1:cidx))]');
                
                %Mat_new = Mat_new + prob.Hes - prob_new.Hes;
                
            end
            
            prob_new.Bk = Mat_new;
            
        end
        
    end
    
    % Original Bfgs update function
%    [ret,prob_new] = bfgsUpdate(newIter,iter,prob_new,prob,par);
    
    % Error
%     Mat = prob_new.Hes + prob_new.Bk;
%     if par.done_correction ~= 0    
%         priReg  = par.last_reg;
%         err     = norm((prob_new.BkC+priReg.*In)-Mat,'fro');
%     else
%         err     = norm(prob_new.BkC-Mat,'fro');
%     end
    
    if ret~=1
        
        ctime = toc(ctime);
        
        break;
    end
    
    g = prob_new.obj_grad;
    
    %%% check residuals
    residual = norm(g,Inf);
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

sol         = iter;
sol.nIter   = nIter;
sol.name    = 'probname';
sol.zoomLS  = par.zoomLS;
sol.obj     = prob.obj;
sol.numFact = par.numFact;
% sol.num_LS = num_LS;
if (bool_conv == 1); sol.ctime = toc(ctime); else sol.ctime = ctime; end
%sol.nhes    = nhes;

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

if storeDeltas == true
    deltas(1,1) = nIter;
    sol.deltas  = deltas;
end

