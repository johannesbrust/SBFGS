function [g] = LIBSVM_grad_VEC(w,y,X)
%------------------------ LIBSVM_grad_VEC --------------------------------%
%
% LIBSVM_grad.m computes the logistic regression gradient value.
%
% INPUTS:
% w := Weights  (size n)
% y := Labels   (size m)
% X := Features (size m x n)
%
% OUTPUTS:
% g := Gradient 
%-------------------------------------------------------------------------%
% 10/04/19, J.B.
% 10/15/19, J.B., vectorization

  % Global variable, regularization parameter.
  global LAM;
  
  [m,n]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
      
  DexpX       = spdiags((-y.*expXwy./(1+expXwy)),0,m,m);
    
  % Gradient computation
  gut        = sum(DexpX*X,1);
    
  g           = gut' + LAM.*w(:,1);
    
end