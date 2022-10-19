function EXPERIMENTS_II_B_FUNC(saveFiles,m,probPath,dataPath,nsol,selectp)
% EXPERIMENTS_II_B_FUNC: Experiments II make comparisons with full-memory
% S-BFGS with varying memory. Uses S-BFGS-P
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
% saveFiles := Flag for saving
% m         := Memory
% probPath  := Problem list
% dataPath  := Path to store outputs
% nsol      := Number of solvers
% selectp   := Indicies to select problems
%-------------------------------------------------------------------------%
% 11/20/19, J.B.

%% Problems

fid         = fopen(probPath);

% Problem file with 13 columns (1st name, 2-13 numeric data). 
pdata       = textscan(fid,'%s%d%d%d%d%d%d%d%d%d%d%d%d'); 


%load('lsidx.mat');
plist       = pdata{1}(:);        

nsel        = length(selectp);
selplist    = plist(selectp);

% Output containers with outputs for selected problems. 
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
par.sbfgsAlg            = 2;
par.addUnknown          = 'setRatio';
par.checkInertia        = 1;%0, 1
par.m                   = m; % 6, 8, 16

par.PrintLV             = 0;
par.MaxIter             = 1000;
par.withHess            = 1;
par.exactRatio          = 0.5;

% BFGS 
parSBFGS                 = get_options();
parSBFGS.alg             = 's-bfgs'; % 's-bfgs', 'exact'
parSBFGS.sbfgsAlg        = 2;
parSBFGS.addUnknown      = 'setRatio';
parSBFGS.checkInertia    = 1;%0, 1
parSBFGS.PrintLV         = 0;
parSBFGS.MaxIter         = 1000;
parSBFGS.exactRatio      = 0.5;

parSBFGS.withHess        = 1;

%% Loop and output store
for j = 1:nsel

    file                   = selplist{j};
    
    prob.name                   = file;
    [~,prob,~,rawfunc]    = read_cutest_prob(file,prob,par);

    clc;
    fprintf('Running Problem: %s \n',file);
    
    n           = prob.n;
    par.spIn    = speye(n);    
    prob.Pmat   = par.spIn;
            
    func        = user_func_setRatio_EXTND(rawfunc);
        
    for k = 1:nsol
        
        if k <= 4
            
            par.Scaling   = k; 
            outData{j,k}  = EX_CSBP_SCALE_F(prob,func,par);
            
        else
            
            outData{j,k}  = EX_SOLVENLP_F(prob,func,parSBFGS);
            
        end
        
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
    
    dname = [dataPath,'experiments','_II_B','_m',num2str(m)];
    save(dname,'outData','outIts','outObjs','outNgs','outTimes');

end

