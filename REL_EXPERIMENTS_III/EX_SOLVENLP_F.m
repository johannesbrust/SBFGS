function [ outexf ] = EX_SOLVENLP_F(prob, fncs, opts)
%---- EX_SOLVENLP_F: Experiment of Compact-Structured L-BFGS Minus V3 func.%
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
% 10/21/19, J.B., solvenlp call

    n   = size(prob.x0,1);
    
    fprintf(' Running: solvenlp on    n=%i \n',n);
        
    [bool_conv,sol] = solvenlp(prob,opts,fncs);
            
    ng          = norm(fncs.getGrad(sol,prob,opts,fncs.rawfunc),Inf);
    
    outexf      = sol;
    outexf.ng   = ng;
    outexf.OK   = bool_conv;
    
    fprintf('Finished: solvenlp in time=%4.2f \n',sol.ctime);
    
    return;

end

