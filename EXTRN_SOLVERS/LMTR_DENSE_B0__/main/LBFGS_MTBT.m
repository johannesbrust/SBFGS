function [x,f,outinfo] = LBFGS_MTBT(fun,x0,params)
% LBFGS_MTBT - limited memory line-search algorithm L-BFGS based on the
%              More-Thuente line search and the initial step is obtained 
%              using backtrack.
%
% For details, see: 
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Linköping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359
% 
% [x,f,outinfo] = LBFGS_MTBT(fun,x0,params) finds a minimizer x of function 
% defined by a handle fun starting from x0 using L-BFGS. f is the function 
% value defined at x.
% The initial parameters are defined in a struct params. 
% The output information is contained in a struct outinfo. 
%
% params contains the following fields {default values}:
%   m       - maximum number of stored vector pairs (s,y) {5}
%   gtol    - tolerance on L2-norm of gradient ||g||<gtol*max(1,||x||) {1e-5}
%   maxit   - maximum number of iterartions {100000}
%   dflag   - display parameter {1}:
%               1 - display every iteration;
%               0 - no display.
%
% outinfo contains the following fields:    
%   ex      - exitflag:
%               1 - norm of gradient is too small
%              -1 - step length is too small
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

% Set parameters for More-Thuente line search procedure cvsrch
gtolls=0.9;
ftolls=1e-4;
xtolls=1e-16;
maxfev=20;
stpmin=1e-20;
stpmax=1e20;

trrad_min = 1e-15; 

% Set linsolve options for solving triangular linear systems
opts1.LT=true;
opts1.RECT=false;
opts1.TRANSA=true;

opts2=opts1;
opts2.TRANSA=false;

% Start measuring CPU time
t0=tic; 

%% Memory Allocation 

n=size(x0,1);

% Allocate memory for elements of L-BFGS matrix B = delta*I + V*W*V', V=[S Y]
S=zeros(n,m);
Y=zeros(n,m);
Sg=zeros(m,1);  % product S'*g for the new iterate
Sg_old=zeros(m,1);  % product S'*g for the previous iterate
Yg=zeros(m,1);  % product Y'*g for the new iterate
Yg_old=zeros(m,1);  % product Y'*g for the previous iterate
E=zeros(m,1);  % E=diag(S'*Y)
L=zeros(m,m);  % transposed upper-triangular part of S'*Y
YY=zeros(m,m);  % Gram matrix Y'*Y

% Initialize indexes and counters
numsy=0;  % number of stored couples (s,y)
maskV=[];  % V(:,maskV)=[S Y]
tract=0;  % number of iterations when TR is active
trrej=0;  % number of iterations when initial trial step is rejected
it=0;  % number of TR iterations
numf=0;  % number of function evaluations
numg=0;  % number of gradient evaluations

%% Initial check for optimality

%   Evaluate function in the starting point
[f0, g0]=fun(x0);
numf=numf+1;
numg=numg+1;
ng=norm(g0);         % L2-norm of gradient

if dflag==1
    fprintf('\n**********************\nRunning L-BFGS-MTBT\n');    
    fprintf('it\t obj\t\t norm(df)\t norm(dx)\t step\n');
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

if ns~=1
    tract=tract+1;
end

% Display information about the last iteration
if (dflag==1)
    fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e\n',it,f,ng,ns,ns);
end 

%% Main loop 
while (mflag==1) 
    %% L-BFGS Update    
    if (sy>1.0e-8*ns*ny)  % positive curvture condition holds, update
        gamma=sy/ny^2;
        if numsy<m  % keep old pairs and add from the new iterate             
            numsy=numsy+1;
            E(numsy)=sy;                          
            maskV=1:numsy;
            indYg=1:numsy-1;
        else  % remove the oldest pair
            E=[E(2:m); sy];            
            YY(1:m-1,1:m-1)=YY(2:m,2:m);
            L(1:m-1,1:m-1)=L(2:m,2:m);
            maskV=maskV([ 2:m 1]);
            indYg=2:m;
        end        
        S(:,maskV(numsy))=s;
        Y(:,maskV(numsy))=y;
        Yg(1:numsy)=Y(:,maskV)'*g; 
        Sg(1:numsy)=S(:,maskV)'*g;
        YY(1:numsy-1,numsy)=Yg(1:numsy-1)-Yg_old(indYg);        
        YY(numsy,1:numsy-1)=YY(1:numsy-1,numsy)';       
        YY(numsy,numsy)=ny^2;
        L(numsy,1:numsy-1)=Sg(1:numsy-1)-Sg_old(indYg);        
        L(numsy,numsy)=sy;
        Yg_old(1:numsy) = Yg(1:numsy); 
        Sg_old(1:numsy) = Sg(1:numsy); 
    else  % skip L-BFGS update but compute S'*g
        Sg_old(1:numsy)=S(:,maskV)'*g; 
    end  % L-BFGS update
    
    %% Quasi-Newton direction
    % Compute the quasi-Newton direction using the inverse Hessian
    % representation by Byrd, Nocedal and Schnabel, 1994    
    w=linsolve(L(1:numsy,1:numsy),Sg_old(1:numsy),opts1);
    w2=linsolve(L(1:numsy,1:numsy), E(1:numsy).*w+...
        gamma*( YY(1:numsy,1:numsy)*w-Yg_old(1:numsy) ), opts2);
    dN=-gamma*g+Y(:,maskV)*(gamma*w) - S(:,maskV)*w2;
    
    % Line search along the quasi-Newton direction    
    g_old=g;
    [x,f,g,step,exls,numfls] = ...
        cvsrch(fun,n,x,f,g_old,dN,1,ftolls,gtolls,xtolls,stpmin,stpmax,maxfev);
    numf = numf + numfls;
    numg = numg + numfls;
    
    if step~=1;
        tract=tract+1;
    end
    
    if exls>1  % line search failed
        ex=exls;
        tcpu=toc(t0);
        outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params);
        return;
    end  
    
    % Compute the trial step    
    s=step*dN;
    ns=norm(dN)*step;
    ng=norm(g);        
    it=it+1;
    
    % Display information about the last iteration
    if dflag==1
        fprintf('%d\t %.3e\t %.3e\t %.3e\t %.3e\n',...
            it,f,ng,ns,step);                
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
