function [varargout] = QUAD_fun( x, g, Q1, Q2 )
%--------------------------- QUAD_fun ------------------------------------%
%
% QUAD_fun.m is a function to compute the structured quadratic objectives:
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
% f     := Function
% gout  := Gradient (if desired)
%-------------------------------------------------------------------------%

  Q   = Q1 + Q2;    
  Qx  = Q*x;
  f   = g'*x + 0.5* (x'*Qx);
    
  varargout{1} = f; 
  
  if nargout == 2
    
    varargout{2} = g + Qx;
    
  end