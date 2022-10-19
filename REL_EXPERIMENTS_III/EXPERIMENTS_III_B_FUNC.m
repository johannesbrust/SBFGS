function EXPERIMENTS_III_B_FUNC(saveFiles,dataPath,N,xl,xu,yl,yu,a,b,sig,nsol)
% EXPERIMENTS_III_B_FUNC: S-BFGS solvers on PDE data.
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
% INPUTS:
% saveFiles := Flag to save files
% dataPath  := Path for storing
% N         := Grid size
% xl        := lower x (2D grid)
% xu        := upper x (2D grid)
% yl        := lower y (2D grid)
% yu        := upper y (2D grid)
% a         := center x direction
% b         := center y direction
% sig       := standard deviation/spread
% nsol      := Number of solvers
%-------------------------------------------------------------------------%
% Initial version: J.B., 11/20/19


% N   = [20,30,40,50,60,70,80,90,100];
% xl  = 0;
% xu  = 1;
% yl  = 0;
% yu  = 1;
% a   = 0.0;
% b   = 0.0;
% sig = 1.0; % 1.0, 5.0

lNs = length(N);
plist = cell(lNs,1);

% Generate multiple problems
for ii=1:lNs

    Nn = N(ii);
    [ yh, g, A] = generate_PE( Nn, xl, xu, yl, yu, a, b, sig );

    p1.yh   = yh;
    p1.g    = g;
    p1.A    = A;
    p1.h    = 1/(Nn*Nn);
    p1.Nn   = Nn;

    p1.name ='Poisson';

    plist{ii,1}   = p1;        

end    
%np = length(plist);

selectp     = 1:length(plist);
nsel        = length(selectp);
selplist    = plist(selectp);

% Output containers with outputs for selected problems. 
%nsol        = 5;
outData     = cell(nsel,nsol);
probData    = cell(nsel,1);

outIts      = zeros(nsel,nsol);
outObjs     = zeros(nsel,nsol);
outNgs      = zeros(nsel,nsol);
outTimes    = zeros(nsel,nsol);
% outStats    = zeros(nsel,nsol);
% outNSkip    = zeros(nsel,nsol);

%probData    = zeros(nsel,4);

%% Options for solvers

% LMTR

% params=struct;
% params.m = 8;  % number of L-BFGS updates
% params.gtol = 1e-6;  % exit if ||g||_2<gtol*max(1,||x||_2)
% params.ranktol = 1e-7;  % tolerance for establishing rank of V
% params.dflag = 0;  % display parameter, 1 if to display information per iteration
% params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
% params.ftol=1e-11;  % tolerance on relative function reduction
% 
% params.maxit = 10000;
% params.thetaT = 1e-3;
% 
% global LAM;
% LAM = 1e-3; % 1e-4 10000

% BFGS 
% parBFGS                 = get_options();
% parBFGS.alg             = 'bfgs'; % 's-bfgs', 'exact'
% parBFGS.sbfgsAlg        = 1;
% parBFGS.addUnknown      = 'setRatio';
% parBFGS.checkInertia    = 0;%0, 1
% parBFGS.PrintLV         = 0;
% parBFGS.MaxIter         = 10000;
% parBFGS.exactRatio      = 0;
% 
% % SBFGS 
% parSBFGS                = get_options();
% parSBFGS.alg            = 's-bfgs'; % 's-bfgs', 'exact'
% parSBFGS.sbfgsAlg       = 1;
% parSBFGS.addUnknown     = 'setRatio';
% parSBFGS.checkInertia   = 0;%0, 1
% parSBFGS.PrintLV        = 0;
% parSBFGS.MaxIter        = 10000;
% parSBFGS.exactRatio     = 0;
% parSBFGS.withHess       = 1;

% CSBMV3 (Compact limited-memory)
par                     = get_options();
par.alg                 = 's-bfgs'; % 's-bfgs', 'exact'
par.sbfgsAlg            = 1;
par.addUnknown          = 'setRatio';
par.checkInertia        = 0;%0, 1
par.m                   = 8; % 6, 8, 16

par.PrintLV             = 0;
par.MaxIter             = 1000;
par.withHess            = 0;

% Regularization
par.alpha               = 0.1; % 0.1

% BFGS 
parBFGS                 = get_options();
parBFGS.alg             = 'bfgs'; % 's-bfgs', 'exact'
parBFGS.sbfgsAlg        = 1;
parBFGS.addUnknown      = 'setRatio';
parBFGS.checkInertia    = 0;%0, 1
parBFGS.PrintLV         = 0;
parBFGS.MaxIter         = 1000;
parBFGS.exactRatio      = 0;

parBFGS.alpha           = par.alpha;

%% Loop and output store
for j = 1:nsel

    file                   = selplist{j};
    prob.name              = file.name;
    
    probData{j}            = file;
    % Vectorization
%     fncLMTR.obj         = @(x)( LIBSVM_fun_VEC(x,y,X));
%     fncLMTR.grad        = @(x)( LIBSVM_grad_VEC(x,y,X));  
%     
%     fncSLMTR.obj        = @(x)( LIBSVM_fun_VEC(x,y,X));
%     fncSLMTR.grad       = @(x)( LIBSVM_grad_STRU_VEC(x,y,X));  
%     fncSLMTR.hprod      = @(x,s)( LIBSVM_hprod_STRU(s));  

    clc;
    fprintf('Running Problem: %s \n',prob.name);
    
    n            = size(file.yh,1);
    prob.x0      = ones(n,1);
    
    prob.n       = n;
    par.spIn     = speye(n);    
    prob.Pmat    = par.spIn;    
    fncSB        = user_func_PDE(file.Nn,file.yh,file.g,file.A,prob.Pmat);
    
    par.Scaling  = 1; 
    out1         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 2;
    out2         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 3;
    out3         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    par.Scaling  = 4;
    out4         = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
    
    out5         = EX_SOLVENLP_F(prob,fncSB,parBFGS);
    
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
    
end

if saveFiles == 1
    
    %dname = './data/experiments_III_B';
    dname = [dataPath,'experiments','_III_B'];
    save(dname,'outData','outIts','outObjs','outNgs','outTimes');

end

