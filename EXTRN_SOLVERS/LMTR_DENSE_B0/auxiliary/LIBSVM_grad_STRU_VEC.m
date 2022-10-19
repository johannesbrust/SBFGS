function [g,gu] = LIBSVM_grad_STRU_VEC(w,y,X)
%------------------------ LIBSVM_grad_STRU_VEC ---------------------------%
%
% LIBSVM_grad_STRU_VEC.m computes the logistic regression structured 
% gradients.
%
% INPUTS:
% w := Weights  (size n)
% y := Labels   (size m)
% X := Features (size m x n)
%
% OUTPUTS:
% g := Gradient (full)
% gu:= Gradient (unknown Hessian part)
%-------------------------------------------------------------------------%
% 10/04/19, J.B.
% 10/15/19, J.B. Structured gradient

  % Global variable, regularization parameter.
  global LAM;
  
  [m,n]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
  
  DexpX     = spdiags((-y.*expXwy./(1+expXwy)),0,m,m);
    
  % Gradient computation
  gut       = sum(DexpX*X,1);
  
  gu        = gut';
    
  g         = gu + LAM.*w(:,1);
              
end