function [ outexf ] = EX_CSBMSV3_F_MF_SCALE(prob, fncs, opts)
%---- EX_CSBMSV3_F_MF_SCALE: Experiment of Compact-Structured L-BFGS Minus 
% V3 matrix-free with Scaling %
%
% This is used for comparison of the compact structured L-BFGS methods to
% external solvers. The proposed methods are an extension to 
% structured BFGS formulas from Petra, Chiang, Anitescu 2018. 
% The extended methods use compact representations of the
% quasi-Newton formulas. The compact forms are
%
% B = B0 - Psi M Psi',
%
% where typically B0 = gamma.I (n x n) multiple of identity initial matrix,
% Psi (n x 2m), M (2m x 2m) are small low rank updates. 
%
% Inputs:
%   prob: Struct as returned by read_cutest_prob to define problem
%   fncs: Struct to compute objective and gradient values
%   opts: Struct to define optimization options
%
% Outputs:
%   outexf: Struct with output information, includes at least objective 
% function value.
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leyffer.
%-------------------------------------------------------------------------%
%
% Initial version: J.B., 09/05/19
% 10/16/19, J.B., Matrix-free implementation
% 10/22/19, J.B., Scaling method

    n   = size(prob.x0,1);
    
    fprintf(' Running: CSBMSV3_MF_SCALE on    n=%i \n',n);
        
    [bool_conv,sol] = CSBM3_SP_MF_SCALE(prob,opts,fncs);
            
    ng          = norm(fncs.getGrad(sol,prob,opts,fncs.rawfunc),Inf);
    
    outexf      = sol;
    outexf.ng   = ng;
    outexf.OK   = bool_conv;
    
    fprintf('Finished: CSBMSV3_MF_SCALE in time=%4.2f \n',sol.ctime);
    
    return;

end

