%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% call user defined function to get derivatives
% 
function prob = evalFunc(iter,prob_in,par,func,hesflag)

 x           = iter.x;
 prob		 = prob_in;

if  hesflag~=-2
 % compute derivatives
 obj	    				= func.getObj(iter,prob,par,func.rawfunc);
 [obj_grad,obj_gk,obj_gu] 	= func.getGrad(iter,prob,par,func.rawfunc);

 % set current derivatives
 prob.obj			= obj ;
 prob.obj_grad 		= obj_grad ;

 if strcmp(par.alg,'exact')
   prob.obj_grad_u 			= [] ;
 end
 
 if strcmp(par.alg,'bfgs')
   prob.obj_grad_u 			= obj_grad ;
 end

 if strcmp(par.alg,'s-bfgs')
   prob.obj_grad_k 	= obj_gk ;
   prob.obj_grad_u 	= obj_gu ;
 end

 if ~strcmp(par.alg,'bfgs') && hesflag>1
   Hes 				= func.getHess(iter,prob,par,func.rawfunc);
   prob.Hes 		= Hes ;
 else 
   prob.Hes			= -1;
 end
elseif hesflag==-2
 if ~strcmp(par.alg,'bfgs')
   Hes 				= func.getHess(iter,prob,par,func.rawfunc);
   prob.Hes 		= Hes ;
 else 
   prob.Hes			= -1;
 end
end

