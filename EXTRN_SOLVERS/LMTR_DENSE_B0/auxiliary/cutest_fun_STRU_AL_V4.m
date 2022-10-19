function [varargout] = cutest_fun_STRU_AL_V4( x, varargin )
%--- Unconstrained structured augmented Lagrangian objective V4 ----------%
%
% cutest_fun_STRU_AL_V4.m is a function to compute the Augmented Lagrangian 
% to interface with the structured SLMTR_DENSE_B0 solver. It implements
% the structured secant approximation from Dennis et al., '89.
%
% L = f(x) + 0.5MU_AL* c(x)^T*c(x),
%
% where f is the CUTEst objective (known Hessian) and c(x)^T*c(x) has
% an unknown Hessian. This function includes a global
% variable MU_AL. 
%
% This function mixes the notion of the 'unknown gradient' (grad(c'c)) 
% and the entire gradient of L according to the paper by Dennis. Thus the 
% outputs are 'arbitrarily' gu, gv4
%
% Inputs: Up to two inputs
% x: Current iterate
% varargin:     - [] empty,
%               - 'Gradient' (compute gradients)
%
% Outputs: Up to three outputs
% varargout:    - Augmented Lagrangian function value,
%               - Known gradient
%               - gu
%               - gv4
%   
%-------------------------------------------------------------------------%
% 09/10/19, J.B., initial version.
% 09/30/19, J.B., modification for Dennis' formula.

  if nargout > 4
      error( 'obj: too many output arguments' );
  end

  % Global variable, penalty parameter.
  global MU_AL;
  
  % Evaluate constraint. Needed for augmented lagrangian objective,
  % and gradient.
  c   = cutest_cons(x);  
  
  if nargin == 1
      
      %s     = varargin{2};
      
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
        
      end
                 
%       else
%         
%         varargout{1}    = ff + 0.5*MU_AL*(c'*c);
%         
%         [~,gk]          = cutest_obj(x);
%         gu              = MU_AL.*cutest_jtprod(x,c);
%         varargout{2}    = gk +  gu;
%         varargout{3}    = gu;
%         
%       end
    % Check if step is included as input
  elseif nargin == 2
      
      ff            = cutest_obj(x);
      
      s             = varargin{1}; 
      
      [~,gk]        = cutest_obj(x);            
      gu            = MU_AL.*cutest_jtprod(x,c);
      
      varargout{1}  = ff + 0.5*MU_AL*(c'*c);
      
      varargout{2}  = gk + gu;
      varargout{3}  = gu;
      varargout{4}  = MU_AL.*cutest_jtprod((x-s),c);
      
  else
      % Only gradients
      
      s             = varargin{1}; 
      
      [~,gk]        = cutest_obj(x);            
      gu            = MU_AL.*cutest_jtprod(x,c);
                  
      varargout{1}  = gk + gu;
      varargout{2}  = gu;
      varargout{3}  = MU_AL.*cutest_jtprod((x-s),c);
      
  end
      
  end