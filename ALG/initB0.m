%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% initialize bfgs update matrix
% 
function prob = initB0(prob,par)

 dim = prob.n;
 
 if strcmp(par.addUnknown,'addNonLinObj_Proj')
    dim = par.num_unknown;
 end

 prob.Bk = par.initB0 * eye(dim);

 prob.last_obj_grad_u = 0;
 prob.last_obj_grad_k = 0;



