function [x,f,outinfo] = LMTR_DDOGL(fun,x0,params)
% LMTR_DDL - limited memory trust-region algorithm D-DOGL.
%
% For details, see: 
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Linköping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359
%
% [x,f,outinfo] = LMTR_DDOGL(fun,x0,params) finds a minimizer x of function 
% defined by a handle fun starting from x0 using D-DOGL. f is the function 
% value defined at x.
% The initial parameters are defined in a struct params. 
% The output information is contained in a struct outinfo. 
%
% params contains the following fields {default values}:
%   m       - maximum number of stored vector pairs (s,y) {5}
%   gtol    - tolerance on L2-norm of gradient ||g||_2<gtol*max(1,||x||_2) {1e-5}
%   ftol    - tolerance on relative function reduction {1e-11}
%   maxit   - maximum number of iterartions {100000}
%   ranktol - tolerance for establishing rank of V {1e-5}
%   dflag   - display parameter {1}:
%               1 - display every iteration;
%               0 - no display.
%
% outinfo contains the following fields:    
%   ex      - exitflag:
%               1 - norm of gradient is too small
%              -1 - TR radius is too small
%              -2 - exceeds maximum number of iterations
%              >1 - line search failed, see cvsrch.m
%   numit   - number of succesful TR iterations
%   numf    - number of function evaluations
%   numg    - number of gradient evaluations
%   tcpu    - CPU time of algorithm execution
%   tract   - number of iterations when TR is active
%   trrej   - number of iterations when initial trial step is rejected
%   params  - input paramaters
% 
% See also LMTR_out, svsrch
%
% Last modified - December 14, 2015
%
% This code is distributed under the terms of the GNU General Public
% License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire 
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose.

% Read input parameters
if nargin<3
    params = struct;
end;
  
if isfield(params,'m') && ~isempty(params.m)
    m = params.m;
else
    m = 5;
end;

if isfield(params,'gtol') && ~isempty(params.gtol)
    gtol = params.gtol;
else
    gtol = 1e-5;
end;

if isfield(params,'ftol') && ~isempty(params.ftol)
    ftol = params.ftol;
else
    ftol = 1e-11;
end;

if isfield(params,'maxit') && ~isempty(params.maxit)
    maxit = params.maxit;
else
    maxit = 100000;
end;

if isfield(params,'dflag') && ~isempty(params.dflag)
    dflag = params.dflag;
else
    dflag = 1;
end;

% Set trust region parameters using the following pseudo-code
%
% if ratio>=tau0 
%   trial step is accepted 
% end
% if ratio>=tau1 and ||s*||>=c3*trrad
%   trrad=c4*trrad 
% elseif ratio<tau2
%   trrad=max(c1*||s||,c2*trrad)
% end
% if trrad<trrad_min
%   return
% end
%
% Accept the trial step also if (ft-f)/abs(f)<ftol && (ratio>0)

tau0=0;                     
tau1=0.25;                  
tau2=0.75;                  
c1=0.5;                     
c2=0.25;                    
c3=0.8;
c4=2;
trrad_min = 1e-15; 

% Set parameters for More-Thuente line search procedure cvsrch
gtolls=0.9;
ftolls=1e-4;
xtolls=1e-16;
maxfev=20;
stpmin=1e-20;
stpmax=1e20;

% Set linsolve options for solving triangular linear systems
opts1.UT=true;
opts1.RECT=false;
opts1.TRANSA=false;

opts2=opts1;
opts2.TRANSA=true;

% Start measuring CPU time
t0=tic;     

%% Memory Allocation 

n=size(x0,1);
m2=m*2;

% Allocate memory for elements of L-BFGS matrix B = delta*I + V*W*V'
V=zeros(n,m2);  % V=[S Y]
Vg=zeros(m2,1);  % product V'*g for the new iterate
Vg_old=zeros(m2,1);  % product V'*g for the previous iterate
VV=zeros(m2,m2);  % Gram matrix V'*V
T=zeros(m,m);  % upper-triangular part of S'*Y
L=zeros(m,m);  % lower-triangular part of S'*Y with 0 main diagonal
E=zeros(m,1);  % E=diag(S'*Y)

% Allocate memory for the solution to TR subproblem s*=-alpha*g-V*p
alpha=0;
beta=0;
p=zeros(m2,1);
p0=zeros(m,1);
TiVg=zeros(m,1); % TiVg=inv(T)*Vg

% Initialize indexes and counters
numsy=0;  % number of stored couples (s,y)
numsy2=0;  % numsy2=numsy*2
maskV=[];  % V(:,maskV)=[S Y]
tract=0;  % number of iterations when TR is active
trrej=0;  % number of iterations when initial trial step is rejected
it=0;  % number of TR iterations
numf=0;  % number of function evaluations
numg=0;  % number of gradient evaluations

%% Initial check for optimality

% Evaluate function and gradient in the starting point
[f0, g0]=fun(x0);
numf=numf+1;
numg=numg+1;
ng=norm(g0);

if dflag==1
    fprintf('\n**********************\nRunning D-DOGL\n');
    fprintf('it\t obj\t\t norm(df)\t norm(dx)\t trrad\n');    
    fprintf('%d\t %.3e\t %.3e\t ---\t\t ---\n',0,f0,ng);
end

if (ng<max(1,norm(x0))*gtol)
    ex=1;
    x=x0;
    f=f0;
    tcpu=toc(t0);
    outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
    return;
end

%% Initialization: line search along normalized steepest descent direction
it=1;
x=x0;
f=f0;
g=g0;
d=-g/ng;  
ns=1;
xt=x+ns*d;    
ft=fun(xt);
numf=numf+1;

% Backtracking
if ft<f  % doubling step length while improvement    
    f=ft;
    ns=ns*2;    
    xt=x+ns*d;    
    ft=fun(xt);
    numf=numf+1;
    while ft<f
        f=ft;
        ns=ns*2;
        xt=x+ns*d;
        ft=fun(xt);
        numf=numf+1;
    end
    ns=ns/2;    
else  % halving step length until improvement        
    while ft>=f
        ns=ns/2;        
        if ns<trrad_min
            tcpu=toc(t0);
            ex=-1;
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst);
            return;
        end
        xt=x+ns*d;
        ft=fun(xt);
        numf=numf+1;
    end
    f=ft;       
end  % line search

g_old=g;
s=ns*d;
x=x+s;
g=fun(x,'gradient');
numg=numg+1;
ng=norm(g);

if ng>gtol*max(1,norm(x))  % norm of gradient is too large
    mflag = 1;  % new iteration is to be made
    y=g-g_old;
    ny=norm(y);
    sy=s'*y;    
    
    % Try More-Thuente line search if positive curvature condition fails
    if (sy<=1.0e-8*ns*ny)
        [x,f,g,ns,exls,numfls] = ...
            cvsrch(fun,n,x0,f0,g_old,d,1,ftolls,gtolls,xtolls,stpmin,stpmax,maxfev);
        numf = numf + numfls;
        numg = numg + numfls;        
        
        if (exls>1)  % line search failed
            ex = exls;
            tcpu=toc(t0); 
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst);
            return;
        end;            

        s=x-x0;
        y=g-g_old;
        ny=norm(y);
        sy=s'*y;
    end 
else
    mflag = 0;  % problem is solved
end

trrad=2*ns;  % take TR radius as the doubled initial step length
if ns~=1
    tract=tract+1;
end

% Display information about the last iteration
if (dflag==1)
    fprintf('D-DOGL\t %d\t %.3e\t %.3e\t %.3e\t ---\n',it,f,ng,ns);
end  


%% Main loop 
while (mflag==1)  % main stopping criterion: ng>gtol*max(1,norm(x))     
    %% L-BFGS update    
    if (sy>1.0e-8*ns*ny)  % positive curvature condition holds, update        
        delta=ny^2/sy;           
        if numsy<m  % keep old pairs and add from the new iterate                                  
            maskVg=1:numsy2;
            maskVV=[1:numsy,m+1:m+numsy];
            if numsy>0
                if Nflag==1  % quasi-Newton step was accepted
                    Vs=(-alpha*Vg_old(maskVg)-VV(maskVV,maskVV)*p(1:numsy2));
                else
                    Vs=-alpha*Vg_old(maskVg)-beta*VV(maskVV,maskVV)*p(1:numsy2);         
                end
            end   
            numsy=numsy+1;
            numsy2=numsy*2;
            maskV=[1:numsy, m+1:m+numsy];
        else  % remove the oldest pair            
            maskV=maskV([ 2:m 1 m+2:m2 m+1]);
            maskVg=[2:m m+2:m2];
            if Nflag==1
                Vs=(-alpha*Vg_old(maskVg)-VV(maskVg,:)*p(1:numsy2));
            else
                Vs=-alpha*Vg_old(maskVg)-beta*VV(maskVg,:)*p(1:numsy2);
            end
            E(1:m-1)=E(2:m);
            VV([1:m-1,m+1:m2-1],[1:m-1,m+1:m2-1])=VV([2:m,m+2:m2],[2:m,m+2:m2]);
            T(1:m-1,1:m-1)=T(2:m,2:m);
            L(1:m-1,1:m-1)=L(2:m,2:m);
        end         
        E(numsy)=sy;            
        V(:,maskV(numsy))=s;
        V(:,maskV(numsy2))=y; 
        VV([numsy,m+numsy],numsy)=[ns^2; sy]; 
        if numsy>1
            VV([1:numsy-1 m+1:m+numsy-1],numsy)=Vs;
        end        
        VV([numsy,m+numsy],m+numsy)=[sy;ny^2];        
        Vg(1:numsy2)=V(:,maskV)'*g;
        VV([1:numsy-1 m+1:m+numsy-1],m+numsy)=Vg([1:numsy-1,numsy+1:numsy2-1])-Vg_old(maskVg);
        VV(numsy,[1:numsy-1 m+1:m+numsy-1])=VV([1:numsy-1 m+1:m+numsy-1],numsy);
        VV(m+numsy,[1:numsy-1 m+1:m+numsy-1])=VV([1:numsy-1 m+1:m+numsy-1],m+numsy);
        T(1:numsy,numsy)=VV(1:numsy,m+numsy); 
        L(numsy,1:numsy-1)=VV(numsy,m+1:m+numsy-1); 
        Vg_old(1:numsy2) = Vg(1:numsy2); 
    else  % skip L-BFGS update but compute V'*g        
        Vg(1:numsy2)=V(:,maskV)'*g;
        Vg_old(1:numsy2) = Vg(1:numsy2);
    end  % L-BFGS update
    
    %% Quasi-Newton step
    % Compute the L2 norm of the quasi-Newton step using the inverse Hessian
    % representation by Byrd, Nocedal and Schnabel, 1994    
    % s=-1/delta*(g+V*p), where p=M*V'*g
    
    alpha=1/delta;
    % Calculate inv(T)*Vg for computing p    
    TiVg(1:numsy)=linsolve(T(1:numsy,1:numsy),Vg(1:numsy),opts1);    
    p0(1:numsy)=E(1:numsy).*TiVg(1:numsy)+...
        alpha*(VV(m+1:m+numsy,m+1:m+numsy)*TiVg(1:numsy)-Vg(numsy+1:numsy2));
    p(1:numsy2)=[linsolve(T(1:numsy,1:numsy),p0(1:numsy),opts2);-alpha*TiVg(1:numsy)];     
    gVp=p(1:numsy2)'*Vg(1:numsy2);    
    nst=sqrt(alpha^2*ng^2+2*alpha*gVp+...
        p(1:numsy2)'*(VV([1:numsy,m+1:m+numsy],[1:numsy,m+1:m+numsy])*p(1:numsy2)));    
    
    ratio=0;    
    trrad_old=trrad;
    af=max(1,abs(f)); 
    if nst<=trrad   % quasi-Newton step is inside TR
        
        % Compute the trial point and evaluate function
        st=-alpha*g-V(:,maskV)*p(1:numsy2);        
        xt=x+st;
        [ft, gt]=fun(xt);
        numf=numf+1;
        
        % Compare actual (ared) and predicted (pred) reductions 
        ared=ft-f;
        if abs(ared)/af<ftol
            ratio=1;
        elseif ared<0 
            mred=-0.5*(alpha*ng^2+p(1:numsy2)'*Vg(1:numsy2));
            ratio=ared/mred;                
        else
            ratio=0;
        end 
        
        if ratio>tau0                   
            Nflag=1;  % quasi-Newton step is accepted
        else
            Nflag=0;  % quasi-Newton step is rejected by ratio
        end       
        
        if (ratio<tau1)        
            trrad=min(c1*nst,c2*trrad);
        end  
        
        if trrad<trrad_min  % TR radius is too small, terminate
            ex = -1;
            tcpu=toc(t0);
            outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
            return
        end
    else  % quasi-Newton step is rejected by TR 
        Nflag=0;            
    end  % checking quasi-Newton step
        
    if ~Nflag    
        %% Quasi-Newton step is rejected, compute double-dogleg step
        % Notations: QN - quasi-Newton step, CP - Cauchy point
        idelta=alpha;
        tract=tract+1;
        
        W=[idelta*VV(1:numsy,1:numsy) idelta*L(1:numsy,1:numsy);...
            idelta*L(1:numsy,1:numsy)' -diag(E(1:numsy))];
        trit=0;        
        nsN=nst;  % norm of QN step
        gBg=delta*ng^2-Vg(1:numsy2)'*(W\Vg(1:numsy2));   
        sigma=-ng^2/gBg;  % step along -g
        gBinvg=ng^2*idelta+gVp;
        tau=0.2+0.8*ng^4/(gBinvg*gBg);      
        phi=-tau*sigma*gBinvg-sigma^2*ng^2;                    
        nw2=tau^2*nsN^2+2*tau*sigma*gBinvg+sigma^2*ng^2;
        if (-sigma*ng < trrad) && (nst>trrad_old) 
            sN=-alpha*g-V(:,maskV)*p(1:numsy2);
        else
            sN=st;
        end
        while (ratio<=tau0)             
            if tau*nsN <= trrad  % trial point is the intersection of QN with the TR boundary
                beta=trrad/nsN;
                st=beta*sN; 
                alpha=beta*idelta;                
                alpha2=0;                
            elseif -sigma*ng > trrad  % trial point is CP which is on the TR boundary
                st=(-trrad/ng)*g;
                beta=0;
                alpha=trrad/ng;
                alpha2=-alpha;
            else              
                % QN is outside and CP is inside. Trial point is the 
                % intersection of the dogleg path with TR boundary.                
                theta=(trrad^2-sigma^2*ng^2);
                theta=theta/(phi+sqrt(phi^2+nw2*theta));
                st=(1-theta)*sigma*g+theta*tau*sN; 
                alpha=(theta-1)*sigma+theta*tau*idelta;
                beta=theta*tau;
                alpha2=(1-theta)*sigma;
            end
            nst=trrad;            
            
            % Compute the trial point
            xt=x+st;
            ft=fun(xt);
            numf=numf+1;
            
            % Compare actual (ared) and predicted (pred) reduction
            ared=ft-f;                        
            if abs(ared)/af<ftol
                ratio=1;
            elseif ared<0
                mred=-alpha*ng^2-beta*gVp+0.5*(alpha2^2*gBg+2*alpha2*beta*ng^2+beta^2*gBinvg);
                ratio=ared/mred;                
            end
            
            if (ratio<tau1)        
                trrad=min(c1*nst,c2*trrad);                
            end  
            
            if trrad<trrad_min  % TR radius is too small, terminate
                ex = -1;
                tcpu=toc(t0); 
                outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
                return
            end
            
            trit=trit+1;
        end  % solving the TR subproblem
        
        % Update counter if the initial trial step is rejected and TR
        % subproblem is to be solved again for a reduced TR radius
        if trit>1
            trrej=trrej+1;
        end        
    end  
    %% Trial step is accepted, update TR radius and gradient
    if (ratio>tau2) && (nst>=c3*trrad)   
        trrad=c4*trrad;   
    end                
    s=st;
    ns=norm(s);
    x=xt;
    f=ft;    
    g_old=g;
    if Nflag==1
        g=gt;
    else
        g=fun(x,'gradient');
    end        
    numg=numg+1;
    ng=norm(g);
    it=it+1; 
    
    % Display information about the last iteration
    if dflag==1 
        fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e\n',...
            it,f,ng,nst,trrad_old);                
    end    
    
    % Check the main stopping criteria
    if ng>gtol*max(1,norm(x))   
        y=g-g_old;
        ny=norm(y);
        sy=s'*y;
    else
        mflag=0;
    end
    
    % Check if the algorithm exceeded maximum number of iterations
    if it > maxit           
        ex=-2;
        tcpu = toc(t0);
        outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
        return;
    end   
    
end  % main loop

ex=1;
tcpu=toc(t0); 
outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
