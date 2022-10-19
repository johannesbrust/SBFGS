function [varargout] = LIBSVM_fun( w, y, X )
%------------------------ LIBSVM_fun -------------------------------------%
%
% LIBSVM_fun.m computes the logistic regression function value,
% and gradient (if desired) for global parameter LAM.
%
% f = 0.5LAM \|w\|^2 + sum_i(log(1 + exp(-w'x_i))),
%
% g = LAM.w          + sum_i(-y_i exp(-w'x_i)/(1+exp(-y_i w'x_i)) )x_i,
%
% where X = [x_1,...x_i] is a feature matrix and y_i are corresponding
% labels.
%
% INPUTS:
% w := Weights  (size n)
% y := Labels   (size m)
% X := Features (size m x n)
%
% OUTPUTS:
% f := Function
% g := Gradient (if desired)
%-------------------------------------------------------------------------%

  % Global variable, regularization parameter.
  global LAM;
  
  [m,n]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
  
  f         = 0.5*LAM*(w'*w) + sum(log(1 + expXwy));
  
  varargout{1} = f; 
  
  if nargout == 2
      
    g(:,1) = LAM.*w(:,1);
    
    for i = 1:m
      
        g(:,1) = g(:,1) + (-y(i)*expXwy(i)/(1+expXwy(i))).*X(i,:)'; 
        
    end
    
    varargout{2} = g;
    
  end