%------------------------------ EX_QUAD_COMP -----------------------------%
% EX_QUAD_COMP Testing IPOPT, LBFGSB and Structured BFGS 
% on structured quadratic functions. Intended as experiment for manuscript. 
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
% Initial version: J.B., 09/05/19
% 09/06/19, J.B., Inclusion of all solvers for single problem.
% 11/06/19, J.B., Comparison on quadratics

%% Resolving dependencies

clc
clear

saveFiles = 1;
fname     = 'EX_QUAD_COMP';

format long e
% warning off all
warning('off','backtrace');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../EXTRN_SOLVERS/IPOPT'));
addpath(genpath('../EXTRN_SOLVERS/LBFGSB'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));

%% Problems

n       = 1000;
g_      = randn(n,1); 
r       = 20;

if n <= r
    r = floor(n/2);
end

lamu    = 999;
laml    = 0;

lamu2   = 1*lamu;
laml2   = 1*laml;

% Forming Q matrices
q1o     = q1_orthog_data(n,r,laml,lamu,1);
q2o     = q1_orthog_data(n,r,laml2,lamu2,1);
%q1d     = q1_diag(n,-0.5*(laml+1),maxd);
%Q2_     = randn(n,n);

% Storing data
p1.g    = g_;

%p1.Q1   = Q1_'*Q1_;
scale   = 1/(1); % 1 % /n, /(n*n)

p1.Q1   = 0.5*(q1o+q1o');
p1.Q2   = scale*0.5*(q2o+q2o');
%p1.Q2   = scale.*Q2_'*Q2_;

p1.name ='P1';
plist   = {p1};

%np = length(plist);

selectp     = 1:length(plist); % 1:10;
nsel        = length(selectp);
selplist    = plist(selectp);

% Output containers with outputs for selected problems. 
nsol        = 4;
outData     = cell(nsel,nsol);   
outProbs    = cell(nsel,1);
outIts      = zeros(nsel,nsol);
outObjs     = zeros(nsel,nsol);
outNgs      = zeros(nsel,nsol);
outTimes    = zeros(nsel,nsol);

%% Options for solvers
par                     = get_options();
par.alg                 = 's-bfgs'; % 's-bfgs', 'exact'
par.sbfgsAlg            = 1;
par.addUnknown          = 'setRatio';
par.checkInertia        = 0;%0, 1
par.m                   = 8; % 6, 8, 16

par.Tol = 5*par.Tol;

par.PrintLV             = 2;
par.MaxIter             = 10000;
par.withHess            = 0;

parSBFGS                = par;
parSBFGS.withHess       = 1;

% BFGS 
parBFGS                 = get_options();
parBFGS.alg             = 'bfgs'; % 's-bfgs', 'exact'
parBFGS.sbfgsAlg        = 1;
parBFGS.addUnknown      = 'setRatio';
parBFGS.checkInertia    = 0;%0, 1
parBFGS.PrintLV         = 0;
parBFGS.MaxIter         = 10000;
parBFGS.exactRatio      = 0;

parBFGS.Tol = 5*parBFGS.Tol;

% % SBFGS
% optsSB                  = get_options();
% optsSB.alg              = 's-bfgs'; % 's-bfgs', 'exact'
% optsSB.sbfgsAlg         = 1;
% optsSB.addUnknown       = 'setRatio';
% optsSB.checkInertia     = 0;
% optsSB.m                = 8; 
% optsSB.mu               = 10000;
% optsSB.PrintLV          = 1;
% optsSB.MaxIter          = 200;

% IPOPT
optsIP.ipopt.hessian_approximation      = 'limited-memory';
optsIP.ipopt.tol                        = 5*1e-6;
optsIP.ipopt.max_iter                   = 10000;
optsIP.ipopt.limited_memory_max_history = 8;
optsIP.ipopt.print_level                = 0;

% LBFGSB
optsLB.m              = 8;
optsLB.maxIts         = 10000;
optsLB.maxTotalIts    = 10000000;
optsLB.printEvery     = Inf;
optsLB.pgtol          = 5*1e-6;

%% Loop and output store
for j = 1:nsel

    file         = selplist{j};
    prob.name    = file.name;
    
    outProbs{j,1} = file;
    
    clc;
    fprintf('Running Problem: %s \n',prob.name);
    
    g            = file.g;
    Q1           = file.Q1;
    Q2           = file.Q2;
    
    n            = size(g,1);
    prob.x0      = zeros(n,1);
    
    prob.n       = n;
    par.spIn     = speye(n);    
    prob.Pmat    = par.spIn;
    fncSB        = user_func_QuadStru(file.g,file.Q1,file.Q2);
    
    quadObj      = @(x)(x'*g+0.5*((x'*(Q1+Q2))*x));
    quadGrad     = @(x)(g+(Q1+Q2)*x);
    
    fncIP.objective = quadObj;
    fncIP.gradient  = quadGrad;

    fncLB        = {quadObj, quadGrad};
    
    par.Scaling  = 1; 
    out1         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 2;
    out2         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 3;
    out3         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 4;
    out4         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    
    outData{j,1} = out1;
    outData{j,2} = out2;
    outData{j,3} = out3;
    outData{j,4} = out4;
    
    for k = 1:nsol
        
        outIts(j,k)   = outData{j,k}.nIter;
        outObjs(j,k)  = outData{j,k}.obj;
        outNgs(j,k)   = outData{j,k}.ng;
        outTimes(j,k) = outData{j,k}.ctime;
    
    end
        
end

if saveFiles == 1
    
    save(fname,'outData','outIts','outObjs','outNgs','outTimes', 'outProbs');

end

