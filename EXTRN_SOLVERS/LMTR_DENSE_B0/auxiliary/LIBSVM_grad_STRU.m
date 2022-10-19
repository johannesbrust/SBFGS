function [g,gu] = LIBSVM_grad_STRU(w,y,X)
%------------------------ LIBSVM_grad_STRU -------------------------------%
%
% LIBSVM_grad_STRU.m computes the logistic regression structured gradients 
%
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

  % Global variable, regularization parameter.
  global LAM;
  
  [m,n]     = size(X);
  gu        = zeros(n,1);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
          
  for i = 1:m
      
     gu(:,1) = gu(:,1) + (-y(i)*expXwy(i)/(1+expXwy(i))).*X(i,:)';  
              
  end
  
  g(:,1)     = LAM.*w(:,1) + gu(:,1);
    
end