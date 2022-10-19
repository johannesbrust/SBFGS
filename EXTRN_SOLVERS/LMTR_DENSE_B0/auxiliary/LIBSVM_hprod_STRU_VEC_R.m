function [r] = LIBSVM_hprod_STRU_VEC_R(w,y,X)
%------------------------ LIBSVM_hprod_STRU_VEC_R ------------------------%
%
% LIBSVM_hprod_STRU_VEC_R.m computes the logistic regression structured 
% Hessian product in reversed (R) order, i.e., known and unknown parts are 
% switched.
%
% INPUTS:
% w := Weights  (size n)
% y := Labels   (size m)
% X := Features (size m x n)
%
% OUTPUTS:
% r := Hessian product (size n)
%-------------------------------------------------------------------------%
% 10/04/19, J.B.
% 10/15/19, J.B., Structured gradient
% 10/23/19, J.B., Reversal of known and unknown gradients 
  
  [m,~]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
  
  DexpX     = spdiags((((y.^2).*expXwy)./((1+expXwy).^2)),0,m,m);
    
  % Hessian product computation
  r         = X'*(DexpX*Xw);
  
%   gut       = sum(DexpX*X,1);
%   
%   gu        = LAM.*w(:,1);
%     
%   g         = gu + gut';
              
end