function [g,gu] = LIBSVM_grad_STRU_VEC_R(w,y,X)
%------------------------ LIBSVM_grad_STRU_VEC_R -------------------------%
%
% LIBSVM_grad_STRU_VEC_R.m computes the logistic regression structured 
% gradients in reversed (R) order, i.e., known and unknown parts are 
% switched.
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
% 10/15/19, J.B., Structured gradient
% 10/23/19, J.B., Reversal of known and unknown gradients 

  % Global variable, regularization parameter.
  global LAM;
  
  [m,~]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
  
  DexpX     = spdiags((-y.*expXwy./(1+expXwy)),0,m,m);
    
  % Gradient computation
  gut       = sum(DexpX*X,1);
  
  gu        = LAM.*w(:,1);
    
  g         = gu + gut';
              
end