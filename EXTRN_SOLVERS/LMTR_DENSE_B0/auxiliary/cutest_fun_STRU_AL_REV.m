function [varargout] = cutest_fun_STRU_AL_REV( x, varargin )
%---- unconstrained structured augmented Lagrangian objective ------------%
%
% cutest_fun_STRU_AL_REV.m is a function to compute the Augmented Lagrangian to 
% interface with the structured SLMTR solvers. This function reverses
% the order of which Hessian is known, and which is unknown. The augmented
% Lagrangian
%
% L = f(x) + 0.5MU_AL* c(x)^T*c(x),
%
% where f is the CUTEst objective (unknown Hessian) and c(x)^T*c(x) has
% a known Hessian. This function includes a global
% variable MU_AL. 
%
% This function computes the 'unknown gradient' in addition 
% to the entire gradient of L.
%
% Inputs: Up to two inputs
% x: Current iterate
% varargin:     - [] empty,
%               - 'Gradient' (compute gradients)
%
% Outputs: Up to three outputs
% varargout:    - Augmented Lagrangian function value,
%               - Known gradient
%               - Unknown gradient
%   
%-------------------------------------------------------------------------%
% 09/10/19, J.B., initial version.

  if nargout > 3
      error( 'obj: too many output arguments' );
  end

  % Global variable, penalty parameter.
  global MU_AL;
  
  % Evaluate constraint. Needed for augmented lagrangian objective,
  % and gradient.
  c   = cutest_cons(x);  
  
  if nargin == 1
      
      % Compute objective function value
      ff  = cutest_obj(x);      
      
      if nargout == 1
          
          varargout{1}    =  ff + 0.5*MU_AL*(c'*c);
      
      elseif nargout == 2 
         
        % Gradient and function requested         
        varargout{1}    = ff + 0.5*MU_AL*(c'*c);
        
        [~,gk]          = cutest_obj(x);
        gu              = MU_AL.*cutest_jtprod(x,c);
        varargout{2}    = gk +  gu;        
                 
      else
        
        varargout{1}    = ff + 0.5*MU_AL*(c'*c);
        
        [~,gk]          = cutest_obj(x);
        gu              = MU_AL.*cutest_jtprod(x,c);
        varargout{2}    = gk +  gu;
        varargout{3}    = gk;
        
      end
  else  
      % Only gradients requested
      
      [~,gk]        = cutest_obj(x);            
      gu            = MU_AL.*cutest_jtprod(x,c);
      
      varargout{1}  = gk + gu;
      varargout{2}  = gk;
      
  end