function func = user_func_PDE(n,yh,g,A,P)
% Structured PDE constrained objective -----------------------------------%
%
% user_func_PDE.m is a function to compute the structured PDE
% constrained objectives
%
% f(w) = 0.5\| P(y-hat{y}) \|^2 + 0.5*alpha*\| w \|,
%
% where PDE constraints, Square Domian
%
%           Ay  = w,    (Interior)
%           y   = g,    (Boundary)
%
% and where Ay := Laplace(y), and P selects elements for the error
% computation. For structured BFGS regularization term \| w \| is assumed
% to have known Hessian.
%
% INPUTS:
% n     := Mesh points (problem dimension (n-2)^2 for 2D with boundaries) 
% yh    := Data
% g     := Boundary data ((n-2)^2 x 4), g = [L,T,R,B], i.e. g =
% [Left,Top,..]
% A     := Finite difference matrix
% P     := Selection matrix for errors
%
% OUTPUT: 
% func := Struct for function evaluations.
%
%-------------------------------------------------------------------------%
% Initial version: 11/04/19, J.B.

func.rawfunc    	= [];
func.getObj         = @(iter,prob,par,rawfunc)(getObj(iter,prob,par,rawfunc,n,yh,g,A,P));
func.getGrad    	= @(iter,prob,par,rawfunc)(getGrad(iter,prob,par,rawfunc,n,yh,g,A,P));
func.getKnowGrad    = @(iter,prob,par,rawfunc)(getKnowGrad(iter,prob,par,rawfunc));
func.getUnknowGrad  = @(iter,prob,par,rawfunc)(getUnknowGrad(iter,prob,par,rawfunc,n,yh,g,A,P));
func.getHess        = @(iter,prob,par,rawfunc)(getHess(iter,prob,par,rawfunc));
func.getHessFull    = @(x,par)(getHessFull(x,par));
func.getHessProd    = @(x,s,par)(getHessProd(x,s,par));

func.computeIdx     = @(n)(computeIdx(n));
end

%%%%%% get objective function
function  f = getObj(iter,prob,par,rawfunc,n,yh,g,A,P)

    [IL,IT,IR,IB] = computeIdx(n);
    
    h       = 1/(n*n);    
    
    alpha   = par.alpha;    
    
    w(:,1)  = -h.*iter.x;
    %w       = iter.x;
    
    w(IL,1)   = w(IL,1) + g(:,1); % Left
    w(IT,1)   = w(IT,1) + g(:,2); % Top
    w(IR,1)   = w(IR,1) + g(:,3); % Right
    w(IB,1)   = w(IB,1) + g(:,4); % Bottom
    
%     w(IL)   = -h.*w(IL) + g(:,1); % Left
%     w(IT)   = -h.*w(IT) + g(:,2); % Top
%     w(IR)   = -h.*w(IR) + g(:,3); % Right
%     w(IB)   = -h.*w(IB) + g(:,4); % Bottom
    
    r       = P*(A\w-yh);
    
    f       = 0.5*(norm(r,2)^2 + alpha*norm(iter.x,2)^2);
        
end

%%%%%% get gradient 
function [gout,gk,gu] = getGrad(iter,prob,par,rawfunc,n,yh,g,A,P)

  [IL,IT,IR,IB] = computeIdx(n);
    
  h       = 1/(n*n);    
    
  alpha   = par.alpha;    
  w(:,1)  = -h.*iter.x;
    %w       = iter.x;
    
  w(IL,1)   = w(IL,1) + g(:,1); % Left
  w(IT,1)   = w(IT,1) + g(:,2); % Top
  w(IR,1)   = w(IR,1) + g(:,3); % Right
  w(IB,1)   = w(IB,1) + g(:,4); % Bottom
    
  r       = P*(A\w-yh);
  
  gk      = alpha.*iter.x;
  
  gu      = A\(P'*(-h.*r));
  
  gout(:,1) = gu + gk; 
  
  if strcmp(par.alg,'bfgs') == true
      
      gk = 0;
      gu = gout;
      
  end
    
end    
    
%%%%%% get gradient of known part
function gk = getKnowGrad(iter,prob,par,rawfunc)
    
    %gk      = Q1*iter.x;
    alpha   = par.alpha;    
    w(:,1)  = iter.x;
    
    gk      = alpha.*w(:,1);
    
    if strcmp(par.alg,'bfgs') == true
      
      gk = 0;     
      
    end
    
end

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iter,prob,par,rawfunc,n,yh,g,A,P)
    
    [IL,IT,IR,IB] = computeIdx(n);

    h       = 1/(n*n);    
    w(:,1)  = -h.*iter.x;

    w(IL,1)   = w(IL,1) + g(:,1); % Left
    w(IT,1)   = w(IT,1) + g(:,2); % Top
    w(IR,1)   = w(IR,1) + g(:,3); % Right
    w(IB,1)   = w(IB,1) + g(:,4); % Bottom

    r       = P*(A\w-yh);

    gu      = A\(P'*(-h.*r));

    if strcmp(par.alg,'bfgs') == true

        alpha   = par.alpha;
        
        gu = gu + alpha.*(iter.x);
        %gu = gu + Q1*iter.x;

    end
     
end

% get Hessian ( known part only )
% requires spIn (sparse Identity (n x n)
function H = getHess(iter,prob,par,rawfunc)
    H = 0;
    
    if strcmp(par.alg,'s-bfgs') == true &&...
            par.withHess == true
        
        H = (par.alpha).*(par.spIn);
        
    end
    
end

% getHessFull ( known part only )
function H = getHessFull(x,par)
    H = (par.alpha).*(par.spIn);
end

% getHessProd ( known part only )
function Hs = getHessProd(x,s,par)
    Hs = (par.alpha).*s;
end

%------------------------ Auxilliary Functions ---------------------------%
% J.B., 11/04/19

function [IL,IT,IR,IB] = computeIdx(n)
% computeIdx -------------------------------------------------------------%
% Compute indices for boundary data

    nm2 = n-2;

    IL  = 1:nm2;
    IT  = 1:nm2:(nm2^2 - (nm2-1));
    IR  = (nm2*(nm2-1)+1):(nm2^2);
    IB  = nm2:nm2:(nm2^2);

end




