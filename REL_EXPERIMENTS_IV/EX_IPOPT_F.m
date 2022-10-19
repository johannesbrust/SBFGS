function [ outexf ] = EX_IPOPT_F(prob, fncs, opts)
%------------------ EX_IPOPT_F: Experiment IPOPT function ----------------%
% Function that uses the IPOPT algorithm (Waechter et al.)
% 
% This is used for comparison with compact structured L-BFGS methods, which
% are an extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
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

    n   = size(prob.x0,1);
    
    fprintf(' Running: IPOPT on    n=%i \n',n);
        
    [x, info]       = ipopt(prob.x0,fncs,opts);
                        
    ng              = norm(fncs.gradient(x),Inf);
    ff              = fncs.objective(x);
    iterations      = info.iter;
    tcpu            = info.cpu;

    status          = info.status;

    outexf.x        = x;
    outexf.nIter    = iterations;    
    outexf.obj      = ff;
    outexf.ctime    = tcpu;
    outexf.ng       = ng;
        
    if status == 0 || status == 1
        okv = 1;
    elseif (status == 2 || status == 3 || status == 4 || status == 5 || status == -1)

        okv = 2;
    else
        okv = 0;
    end
    
    outexf.OK       = okv;
    outexf.misc1    = status;
    
    fprintf('Finished: IPOPT in time=%4.2f \n',tcpu);
    
    return;

end

