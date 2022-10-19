%------------------------------ EXPERIMENTS_II_A -------------------------%
% EXPERIMENTS_II_A: Memory comparison m=8,50 and original S-BFGS solver.
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. The compact forms are
%
% B = B0 - Psi M Psi',
%
% where typically B0 = gamma.I (n x n) multiple of identity initial matrix,
% Psi (n x 2m), M (2m x 2m) are small low rank updates. 
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leyffer.
%
%-------------------------------------------------------------------------%
%
% Initial version: J.B., 09/10/19
% 09/30/19, J.B., Implementation using V4
% 10/04/19, J.B., Initial test on LIBSVM problem
% 10/07/19, J.B., Including CSBMV3
% 10/15/19, J.B., List of problems, removal of line-search solver,
%                   inclusion of vectorized objectives
% 10/16/19, J.B., Including Matrix-free implementation
% 10/17/19, J.B., Including SLMTR_V2_INTF1
% 10/21/19, J.B., Quadratic tests, and use of solvenlp
% 10/22/19, J.B., Use functions to form Q test matrices
% 10/22/19, J.B., Implementation of scaling tests
%
% (1): delta = bk'*bk/sk'*bk,
% (2): delta = sk'*bk/sk'*sk,
% (3): delta = yk'*yk/sk'*yk,
% (4): delta = sk'*yk/sk'*sk,
%
% where bk = yk + K(xk1)sk, yk = grad(u(xk1)) - grad(u(xk)),
% and K(x) is the known Hessian, and the Hessian of u(x) is unknown.
% 10/24/19, J.B., Including BFGS solver (line-search)
% 11/19/19, J.B., Manuscript run

%% Resolving dependencies

clc
clear

saveFiles = 1;
fname     = 'EXPERIMENTS_II_A';

format long e
% warning off all
warning('off','backtrace');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

%addpath(genpath('../EXTRN_SOLVERS/LMTR_DENSE_B0'));
%addpath(genpath('../ALG'));
addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));

%% Problems

fid                 = fopen('fullProblemTable.txt');

% Problem file with 13 columns (1st name, 2-13 numeric data). 
pdata               = textscan(fid,'%s%d%d%d%d%d%d%d%d%d%d%d%d'); 


%load('lsidx.mat');
plist               = pdata{1}(:);        

%np = length(plist);

selectp     = 2;
%selectp     = 1:length(plist);
nsel        = length(selectp);
selplist    = plist(selectp);

% Output containers with outputs for selected problems. 
nsol        = 5;
outData     = cell(nsel,nsol);   

outIts      = zeros(nsel,nsol);
outObjs     = zeros(nsel,nsol);
outNgs      = zeros(nsel,nsol);
outTimes    = zeros(nsel,nsol);
%outStats    = zeros(nsel,nsol);
%outNSkip    = zeros(nsel,nsol);

%probData    = zeros(nsel,4);

%% Options for solvers

par                     = get_options();
par.alg                 = 's-bfgs'; % 's-bfgs', 'exact'
par.sbfgsAlg            = 1;
par.addUnknown          = 'setRatio';
par.checkInertia        = 0;%0, 1
par.m                   = 8; % 6, 8, 16

par.PrintLV             = 0;
par.MaxIter             = 10000;
par.withHess            = 0;
par.exactRatio          = 0.5;

% BFGS 
parSBFGS                 = get_options();
parSBFGS.alg             = 's-bfgs'; % 's-bfgs', 'exact'
parSBFGS.sbfgsAlg        = 1;
parSBFGS.addUnknown      = 'setRatio';
parSBFGS.checkInertia    = 0;%0, 1
parSBFGS.PrintLV         = 0;
parSBFGS.MaxIter         = 10000;
parSBFGS.exactRatio      = 0.5;

parSBFGS.withHess        = 0;

%% Loop and output store
for j = 1:nsel

    file                   = selplist{j};
    
    prob.name                   = file;
    [ret,prob,par_U,rawfunc]    = read_cutest_prob(file,prob,par);

    clc;
    fprintf('Running Problem: %s \n',file);
    
    n           = prob.n;
    par.spIn    = speye(n);    
    prob.Pmat   = par.spIn;
            
    func        = user_func_setRatio_EXTND(rawfunc);
    
    par.Scaling  = 1; 
    out1         = EX_CSBMSV3_F_MF_SCALE(prob,func,par);
    par.Scaling  = 2;
    out2         = EX_CSBMSV3_F_MF_SCALE(prob,func,par);
    par.Scaling  = 3;
    out3         = EX_CSBMSV3_F_MF_SCALE(prob,func,par);
    par.Scaling  = 4;
    out4         = EX_CSBMSV3_F_MF_SCALE(prob,func,par);
    
    out5         = EX_SOLVENLP_F(prob,func,parSBFGS);
    
%     outLMS2      = EX_LMTR_F(prob,fncSLMTR,params,9);
    
    outData{j,1} = out1;
    outData{j,2} = out2;
    outData{j,3} = out3;
    outData{j,4} = out4;
    outData{j,5} = out5;
    
    for k = 1:nsol
        
        outIts(j,k)   = outData{j,k}.nIter;
        outObjs(j,k)  = outData{j,k}.obj;
        outNgs(j,k)   = outData{j,k}.ng;
        outTimes(j,k) = outData{j,k}.ctime;
        
        %outStats(j,k) = outData{j,k}.misc1;
        %outNSkip(j,k) = outData{j,k}.misc6;
    
    end
    
    cutest_terminate();
        
    % Delete cutest files
    delete *.f;
    delete *.o;
    delete *.dylib;
    delete *.d;
    delete *.mexmaci64;
    
end

if saveFiles == 1
    
    save(fname,'outData','outIts','outObjs','outNgs','outTimes');

end

