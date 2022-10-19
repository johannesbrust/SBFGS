function [g] = LIBSVM_grad(w,y,X)
%------------------------ LIBSVM_grad ------------------------------------%
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

  % Global variable, regularization parameter.
  global LAM;
  
  [m,n]     = size(X);
  
  Xw        = X*w;
  
  expXwy    = exp(-Xw.*y);
      
  g(:,1) = LAM.*w(:,1);
    
  for i = 1:m
      
     g(:,1) = g(:,1) + (-y(i)*expXwy(i)/(1+expXwy(i))).*X(i,:)'; 
        
  end
    
end