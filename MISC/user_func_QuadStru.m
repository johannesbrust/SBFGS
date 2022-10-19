function func = user_func_QuadStru(g,Q1,Q2)
%--------------------- Structured Quadratic objective --------------------%
%
% user_func_QuadStru.m is a function to compute the structured quadratic objectives:
%
% f = x' g + 0.5* x'(Q1 + Q2)x,
%
% where Q2 is the unknown Hessian.
%
% INPUT:
% g := Gradient 
% Q1 := Known Hessian
% Q2 := Unknown Hessian
%
% OUTPUT: 
% func := Struct for function evaluations.
%
%-------------------------------------------------------------------------%
% Initial version: 08/23/19, J.B.
% 10/07/19, J.B., LIBSVM version
% 10/16/19, J.B., Matrix-free version (MF), vectorization
% 10/21/19, J.B., Structured quadratic implementation, modifications for
% algorithms comparisons.
% 10/23/19, J.B., Including g in the definition of known gradient. 

func.rawfunc    	= [];
func.getObj         = @(iter,prob,par,rawfunc)(getObj(iter,prob,par,rawfunc,g,Q1,Q2));
func.getGrad    	= @(iter,prob,par,rawfunc)(getGrad(iter,prob,par,rawfunc,g,Q1,Q2));
func.getKnowGrad    = @(iter,prob,par,rawfunc)(getKnowGrad(iter,prob,par,rawfunc,g,Q1));
func.getUnknowGrad  = @(iter,prob,par,rawfunc)(getUnknowGrad(iter,prob,par,rawfunc,g,Q1,Q2));
func.getHess        = @(iter,prob,par,rawfunc)(getHess(iter,prob,par,rawfunc,Q1));
func.getHessFull    = @(x,par)(getHessFull(x,par,Q1));
func.getHessProd    = @(x,s,par)(getHessProd(x,s,par,Q1));

end

%%%%%% get objective function
function  f = getObj(iter,prob,par,rawfunc,g,Q1,Q2)

    Q = Q1 + Q2;
    x = iter.x;
    
    f = g'*x + 0.5* ((x'*Q)*x);
        
end

%%%%%% get gradient 
function [gout,gk,gu] = getGrad(iter,prob,par,rawfunc,g,Q1,Q2)

  x     = iter.x;
  
  %gk   = Q1*x;  
  gk    = g + Q1*x;
  
  %gu   = g+ Q2*x;
  gu    = Q2*x;
  
  gout(:,1) = gk + gu; 
  
  if strcmp(par.alg,'bfgs') == true
      
      gk = 0;
      gu = gout;
      
  end
    
end    
    
%%%%%% get gradient of known part
function gk = getKnowGrad(iter,prob,par,rawfunc,g,Q1)
    
    %gk      = Q1*iter.x;
    gk      = g + Q1*iter.x;
    
    if strcmp(par.alg,'bfgs') == true
      
      gk = 0;     
      
    end
    
end

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iter,prob,par,rawfunc,g,Q1,Q2)
    
    gu        = Q2*iter.x;
    %gu        = g + Q2*iter.x;
    
    if strcmp(par.alg,'bfgs') == true
      
      gu = gu + g + Q1*iter.x;     
      %gu = gu + Q1*iter.x;
      
    end
     
end

% get Hessian ( known part only )
% requires spIn (sparse Identity (n x n)
function H = getHess(iter,prob,par,rawfunc,Q1)
    H = 0;
    
    if strcmp(par.alg,'s-bfgs') == true &&...
            par.withHess == true
        
        H = Q1;
        
    end
    
end

% getHessFull ( known part only )
function H = getHessFull(x,par,Q1)
    H = Q1;
end

% getHessProd ( known part only )
function Hs = getHessProd(x,s,par,Q1)
    Hs = Q1*s;
end





