function [ outexf ] = EX_LBFGSB_F(prob, fncs, opts)
%------------------ EX_LBFGSB_F: Experiment LBFGSB function --------------%
% Function that uses the L-BFGS-B algorithm (wrapper S.Becker)
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
    
    fprintf(' Running: L-BFGS-B on    n=%i \n',n);
        
    Ons = ones(n,1);
    l   = -Inf.*Ons;
    u   = +Inf.*Ons;

    opts.x0 = prob.x0;
    opts.l  = l;
    opts.u  = u;

    tlbfgs = tic;

    [xopt,fopt,info]    = lbfgsb(fncs,l,u,opts);

    tlbfgs_o = toc(tlbfgs);

    iterations      = info.iterations;
    totalIterations = info.totalIterations;
    lbfgsmes        = info.lbfgs_message1;
    errhist         = info.err;
    if isfield(info,'taskInteger')==true; taskInteger  = info.taskInteger; else taskInteger = 0; end;    
    
    outexf.x        = xopt;
    outexf.nIter    = iterations;
    outexf.iter     = totalIterations;
    outexf.obj      = fopt;
    outexf.ctime    = tlbfgs_o;
    outexf.ng       = norm(fncs{2}(xopt),Inf);
    outexf.misc1    = taskInteger;
    outexf.misc2    = lbfgsmes;
    outexf.misc3    = errhist;
    
    if taskInteger == 21 || taskInteger == 22
        outexf.OK = 1;
    elseif taskInteger == 31 || taskInteger == 32 || taskInteger == 33
        outexf.OK = 2;
    else
        outexf.OK = 0;
    end
    
    fprintf('Finished: L-BFGS-B in time=%4.2f \n',tlbfgs_o);

end

