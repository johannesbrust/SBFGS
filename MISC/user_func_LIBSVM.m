function func = user_func_LIBSVM(y,X)
%--------------------- LIBSVM structured objective -----------------------%
%
% user_func_LIBSVM.m is a function to compute the structured logistic
% regression:
%
% f = 0.5LAM \| w \|^2 + sum_i(log(1 + exp(-w'x_i))),
%
% g = LAM.w          + sum_i(-y_i exp(-w'x_i)/(1+exp(-y_i w'x_i)) )x_i,
%
% where X = [x_1,...x_i] is a feature matrix and y_i are corresponding
% labels.
%
% Here the regularization term 0.5 LAM. \| w \|^2 is assumed to be the
% known part.
%
%
% INPUT:
% y := Labels (LIBSVM)
% X := Features (LIBSVM)
%
% OUTPUT: 
% func := Struct for function evaluations.
%
%-------------------------------------------------------------------------%
% Initial version: 08/23/19, J.B.
% 10/07/19, J.B., LIBSVM version

func.rawfunc    	= [];
func.getObj         = @(iter,prob,par,rawfunc)(getObj(iter,prob,par,rawfunc,y,X));
func.getGrad    	= @(iter,prob,par,rawfunc)(getGrad(iter,prob,par,rawfunc,y,X));
func.getKnowGrad    = @(iter,prob,par,rawfunc)(getKnowGrad(iter,prob,par,rawfunc));
func.getUnknowGrad  = @(iter,prob,par,rawfunc)(getUnknowGrad(iter,prob,par,rawfunc,y,X));
func.getHess        = @(iter,prob,par,rawfunc)(getHess(iter,prob,par,rawfunc));

end

%%%%%% get objective function
function  f = getObj(iter,prob,par,rawfunc,y,X)

    LAM       = par.LAM;
    
    Xw        = X*iter.x;  
    expXwy    = exp(-Xw.*y);
  
    f         = 0.5*LAM*(iter.x'*iter.x) + sum(log(1 + expXwy));
    
end

%%%%%% get gradient 
function [g,gk,gu] = getGrad(iter,prob,par,rawfunc,y,X)

  LAM       = par.LAM;
  
  [m,n]     = size(X);
  gu        = zeros(n,1);
  
  Xw        = X*iter.x;
  
  expXwy    = exp(-Xw.*y);
          
  for i = 1:m
      
     gu(:,1) = gu(:,1) + (-y(i)*expXwy(i)/(1+expXwy(i))).*X(i,:)';  
              
  end
  
  gk        = LAM*iter.x(:,1);
  
  g(:,1)    = gk + gu(:,1); 
    
end    
    
%%%%%% get gradient of known part
function gk = getKnowGrad(iter,prob,par,rawfunc)
    
    LAM     = par.LAM;   
    gk      = LAM*iter.x(:,1);
    
end

%%%%%% get gradient of unknown part
function gu = getUnknowGrad(iter,prob,par,rawfunc,y,X)

    [m,n]     = size(X);
    gu        = zeros(n,1);
  
    Xw        = X*iter.x;
  
    expXwy    = exp(-Xw.*y);
          
    for i = 1:m
      
        gu(:,1) = gu(:,1) + (-y(i)*expXwy(i)/(1+expXwy(i))).*X(i,:)';  
              
    end
     
end

% get Hessian ( known part only )
% requires spIn (sparse Identity (n x n)
function H = getHess(iter,prob,par,rawfunc)
    H = par.LAM*par.spIn;
end
