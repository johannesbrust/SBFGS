function [gout] = QUAD_grad(x, g, Q1, Q2)
% QUAD_grad Quadratic gradient 
%
% QUAD_grad_STRU.m is a function to compute the structured quadratic gradients:
%
% f = x' g + 0.5* x'(Q1 + Q2)x,
%
% where Q2 is the unknown Hessian.
%
% INPUTS:
% x     := Variables        (size n)
% g     := Linear term      (size n)
% Q1    := Known Hessian    (size n x n)
% Q2    := Unknown Hessian  (size n x n)
%
% OUTPUTS:
% gout  := Gradient 
% gu    := Unknown Gradient
%-------------------------------------------------------------------------%
        
  Q1x = Q1*x;
  Q2x = Q2*x;
  
  gu  = Q2x;
  gout= gu + Q1x + g;
  
  