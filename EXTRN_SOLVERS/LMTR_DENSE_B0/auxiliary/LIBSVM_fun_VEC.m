function [varargout] = LIBSVM_fun_VEC( w, y, X )
%------------------------ LIBSVM_fun_VEC ---------------------------------%
%
% LIBSVM_fun_VEC.m computes the logistic regression function value,
% and gradient (if desired) for global parameter LAM using vectorization.
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
      
    DexpX       = spdiags((-y.*expXwy./(1+expXwy)),0,m,m);
    
    % Gradient computation
    g           = sum(DexpX*X,1);
    
    varargout{2} = g' + LAM.*w(:,1);
    
  end