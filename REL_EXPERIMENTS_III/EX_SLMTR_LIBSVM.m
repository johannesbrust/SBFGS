%------------------------------ EX_SLMTR_LIBSVM --------------------------%
% EX_SLMTR_LIBSVM runs SLMTR, which is a trust-region structured 
% L-BFGS method on LIBSVM test problems. 
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

%% Resolving dependencies

clc
clear

saveFiles = 1;
fname     = 'EX_SLMTR_LIBSVM';

format long e
% warning off all
warning('off','backtrace');
warning('off','MATLAB:singularMatrix');
warning('off','MATLAB:nearlySingularMatrix');

addpath(genpath('../EXTRN_SOLVERS/LMTR_DENSE_B0'));
addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));


LIBSVMPath = '/Users/johannesbrust/Dropbox/codes/LIBSVM/';

addpath(genpath(LIBSVMPath));

%probPath = [LIBSVMPath,'libsvm/'];
probPath = [LIBSVMPath,'data/'];

%% Problems

%plist = {'heart_scale'};
plist = {'gisette/gisette_scale',...
            'colon_cancer/colon-cancer',...
            'leukemia/leu',...
            'real_sim/real-sim',...
            'duke/duke.tr',...
            'rcv1/rcv1_train.binary',...
            'a9a/a9a',...
            'mushrooms/mushrooms',...
            'w8a/w8a',...
            'madelon/madelon'};        

%np = length(plist);

selectp     = 1:length(plist);
nsel        = length(selectp);
selplist    = plist(selectp);

% Output containers with outputs for selected problems. 
nsol        = 4;
outData     = cell(nsel,nsol);   

outIts      = zeros(nsel,nsol);
outObjs     = zeros(nsel,nsol);
outNgs      = zeros(nsel,nsol);
outTimes    = zeros(nsel,nsol);
outStats    = zeros(nsel,nsol);
outNSkip    = zeros(nsel,nsol);

%probData    = zeros(nsel,4);

%% Options for solvers

% LMTR

params=struct;
params.m = 8;  % number of L-BFGS updates
params.gtol = 1e-6;  % exit if ||g||_2<gtol*max(1,||x||_2)
params.ranktol = 1e-7;  % tolerance for establishing rank of V
params.dflag = 0;  % display parameter, 1 if to display information per iteration
params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
params.ftol=1e-11;  % tolerance on relative function reduction

params.maxit = 10000;
params.thetaT = 1e-3;

global LAM;
LAM = 1e-3; % 1e-4 10000

% CSBMV3

par                 = get_options();
par.alg             = 's-bfgs'; % 's-bfgs', 'exact'
par.sbfgsAlg        = 1;
par.addUnknown      = 'setRatio';
par.checkInertia    = 0;%0, 1
par.m               = 8; % 6, 8, 16

par.LAM             = LAM;

par.PrintLV         = 0;

par.MaxIter         = 10000;

%% Loop and output store
for j = 1:nsel

    file                   = selplist{j};
    prob.name              = file;
    
    %% Read LIBSVM problem. 
    [y,X] = libsvmread([probPath,file]);
    
    % Vectorization
    fncLMTR.obj         = @(x)( LIBSVM_fun_VEC(x,y,X));
    fncLMTR.grad        = @(x)( LIBSVM_grad_VEC(x,y,X));  
    
    fncSLMTR.obj        = @(x)( LIBSVM_fun_VEC(x,y,X));
    fncSLMTR.grad       = @(x)( LIBSVM_grad_STRU_VEC(x,y,X));  
    fncSLMTR.hprod      = @(x,s)( LIBSVM_hprod_STRU(s));  
    
    clc;
    fprintf('Running Problem: %s \n',file);
    
    [m,n]        = size(X);
    prob.x0      = zeros(n,1);
    
    prob.n       = n;
    par.spIn     = speye(n);    
    prob.Pmat    = par.spIn;
    fncSB        = user_func_LIBSVM_MF(y,X);
    
    outLM        = EX_LMTR_F(prob,fncLMTR,params,8);
    outLMS       = EX_LMTR_F(prob,fncSLMTR,params,7);
    outCSB       = EX_CSBMSV3_F_MF(prob,fncSB,par);    
    outLMS2      = EX_LMTR_F(prob,fncSLMTR,params,9);
    
    outData{j,1} = outLM;
    outData{j,2} = outLMS;
    outData{j,3} = outCSB;
    outData{j,4} = outLMS2;
    
    for k = 1:nsol
        
        outIts(j,k)   = outData{j,k}.nIter;
        outObjs(j,k)  = outData{j,k}.obj;
        outNgs(j,k)   = outData{j,k}.ng;
        outTimes(j,k) = outData{j,k}.ctime;
        
        %outStats(j,k) = outData{j,k}.misc1;
        %outNSkip(j,k) = outData{j,k}.misc6;
    
    end
    
end

if saveFiles == 1
    
    save(fname,'outData','outIts','outObjs','outNgs','outTimes','outStats');

end

