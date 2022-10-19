%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% update B by (6.19) 
%
% 
function [ret,prob_new] = bfgsUpdate(newiter,iter,prob_new_in,prob,par)

alg				= par.alg;
bfgsAlg			= par.sbfgsAlg;
Pmat            = prob.Pmat;

prob_new        = prob_new_in;

dim_x       	= prob.n;
x_new           = newiter.x;
x_old           = iter.x;
grad_new		= prob_new.obj_grad_u;
grad_old		= prob.obj_grad_u;
Mat			  	= prob.Bk;


ret             = 1;

gamma = 0;
if strcmp(alg,'s-bfgs') && par.skipBFGS==0
  sk 		= prob_new.sk;
  yk 		= prob_new.yk;    
  ak 		= prob_new.ak;
  bk 		= prob_new.bk;
  alpha  	= prob_new.alpha;
  gamma  	= prob_new.gamma;
  if gamma <=0 && (bfgsAlg==1 && ~strcmp(par.addUnknown,'addNonLinObj_Proj') ) % for other cases, we can do regularization
    par.skipBFGS=1;
  end  

elseif strcmp(alg,'bfgs')
  sk = x_new - x_old;
  yk = grad_new - grad_old;
  ak = Mat*sk;
  bk = yk;
  alpha  = sk'*ak;
  gamma  = sk'*bk;
  if gamma <=0 
    par.skipBFGS=1;
  end
end

  % update B 

if (par.stopBFGS==1)
  warning('set bfgs approximat matrix to 0');
  Mat_new = zeros(length(sk));
else

%   if (par.skipBFGS==1) && (  strcmp(alg,'bfgs') || ( strcmp(alg,'s-bfgs') && bfgsAlg==1)   )
  if (par.skipBFGS==1) 
    warning('skip current bfgs update');
    ret=0;
	Mat_new = Mat;
	prob_new.Hes = prob.Hes;
	ret=1;
  else
    Mat_new = Mat - (ak*ak')/alpha + bk*bk'/gamma;

    if bfgsAlg==1 && strcmp(alg,'s-bfgs')
        if strcmp(par.addUnknown,'addNonLinObj_Proj')
            Mat_new = Mat_new + Pmat*prob.Hes*Pmat' - Pmat*prob_new.Hes*Pmat';
        else
            Mat_new = Mat_new + prob.Hes - prob_new.Hes;
        end         
    end
  end
  
end

  % check if new hessian approximation is pd
%  if strcmp(alg,'bfgs')
%    Hes_use = zeros(dim_x,dim_x);
%  elseif bfgsAlg==1 && strcmp(alg,'s-bfgs')
%   	Hes_use   = prob.Hes;
%  elseif bfgsAlg==2 && strcmp(alg,'s-bfgs')
%	Hes_use   = prob_new.Hes;
%  else
%   	error('bfgsAlg should be 1 or 2 within s-bfgs ');
%  end 
%
%   oldAppro  = Mat+Hes_use;
%   TestMat1  = (eye(dim_x)-sk*bk'/gamma);
%   TestMat   = TestMat1*(inv(oldAppro))*TestMat1'+(sk*sk')/gamma;
%   newAppro1 = inv(TestMat);
% 
%   if ~strcmp(alg,'bfgs')
%     newAppro2 = Mat_new + prob_new.Hes;
%   else 
%     newAppro2 = Mat_new;  
%   end
%   
%   checek = 11;
%   if  ~all(eig(oldAppro)>0)
%       warning('Current Hessian approximation is not PD');
%   end  
%   if  ~all(eig(newAppro2)>0) 
%       warning('New Hessian approximation is not PD');
%   end
% 
%   if max(max(abs(newAppro1-newAppro2))) >1e-3 
%      warning('nemerical error! cond(mat1) = %2.2e, cond(mat2) = %2.2e, difference: %3.6e', cond(newAppro1),cond(newAppro2),max(max(abs(newAppro1-newAppro2))) );
%   end 

prob_new.Bk = Mat_new;




