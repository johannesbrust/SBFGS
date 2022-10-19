function EXPERIMENTS_IV_FUNC(saveFiles,PHI,datapath,...
    nsol,nruns,ns,rscale,selectp)
% EXPERIMENTS_IV_FUNC: Function to run comparisons on structured Quadratics
% In additon to 4 L-S-BFGS-M and 4 L-S-BFGS-P external L-BFGS-B and IPOPT
% included.
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
% See change-log and variables in the script EXPERIMENT_I for details
% of input arguments in this function.
%
% INPUTS:
% saveFiles := Flag save/not save
% PHI       := LRG/SM clustering
% datapath  := Data storing path
% nsol      := Number of solvers
% nruns     := Number of simulations
% ns        := Vector of problem dimensions
% rscale    := Parameter to determine clustering
% selectp   := Indices of problems to select
%-------------------------------------------------------------------------%
% 11/20/19, J.B.

%% Resolving dependencies

% Phi parameter to cluster eigenvalues
% -SM : Small eigenvalues clustered
% -LRG: Large eigenvalues clustered
if exist('PHI','var') == false
    PHI = 'SM';
elseif isempty(PHI) == true 
    PHI = 'SM';
end

% format long e
% % warning off all
% warning('off','backtrace');
% warning('off','MATLAB:singularMatrix');
% warning('off','MATLAB:nearlySingularMatrix');

rs          = ns./rscale; %10
nsel        = length(ns);

% Eigenvalue parameters
if strcmp(PHI,'LRG') == true
    
    phi         = 1000; %0
    lamu        = -999; % 999
    laml        = 0; % 0

    lamu2       = -999; % 999
    laml2       = 0; % 0
    
else
    
    phi         = 1; %0
    lamu        = 999;
    laml        = 0;

    lamu2       = 999;
    laml2       = 0;
    
end

nsolm2      = nsol-2;

outData     = cell(nruns,nsel,nsol);   

outIts      = zeros(nsel,nsol);
outObjs     = zeros(nsel,nsol);
outNgs      = zeros(nsel,nsol);
outTimes    = zeros(nsel,nsol);
%outStats    = zeros(nsel,nsol);
%outNSkip    = zeros(nsel,nsol);
%outErrs     = zeros(nsel,2*nsol);

%% Options for solvers

% Minus
par                     = get_options();
par.alg                 = 's-bfgs'; % 's-bfgs', 'exact'
par.sbfgsAlg            = 1;
par.addUnknown          = 'setRatio';
par.checkInertia        = 0;%0, 1
par.m                   = 8; % 6, 8, 16
par.Tol                 = 5*par.Tol;

par.PrintLV             = 0;
par.MaxIter             = 10000;
par.withHess            = 0;
par.storeDeltas         = true;

% Plus
parP                     = get_options();
parP.alg                 = 's-bfgs'; % 's-bfgs', 'exact'
parP.sbfgsAlg            = 2;
parP.addUnknown          = 'setRatio';
parP.checkInertia        = 1;%0, 1
parP.m                   = 8; % 6, 8, 16
parP.Tol                 = 5*parP.Tol;

parP.PrintLV             = 0;
parP.MaxIter             = 10000;
parP.withHess            = 1;
parP.storeDeltas         = true;

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
optsLB.factr          = 1e1; % 1e7

% Original BFGS
% parBFGS                 = get_options();
% parBFGS.alg             = 'bfgs'; % 's-bfgs', 'exact'
% parBFGS.sbfgsAlg        = 1;
% parBFGS.addUnknown      = 'setRatio';
% parBFGS.checkInertia    = 0;%0, 1
% parBFGS.PrintLV         = 0;
% parBFGS.MaxIter         = 10000;
% parBFGS.exactRatio      = 0;

%% 'nruns' repeated computations
for ir = 1:nruns

    nps      = length(ns);
    plist    = cell(nps,1);

    for p = 1:nps

        
        %% Problems
        % Quadratics with 'r' eigenvalues within 1 <= lam <= 1000 
        % and the remaining lam = 1.
        n       = ns(p,1);
        r       = rs(p,1);

        g_      = randn(n,1);
        if n <= r
            r = floor(n/2);
        end

        % Forming Q matrices
        q1o     = q1_orthog_data(n,r,laml,lamu,phi);
        q2o     = q1_orthog_data(n,r,laml2,lamu2,phi);

        %Q2_     = randn(n,n);

        % Storing data
        p1.g    = g_;

        p1.Q1   = 0.5*(q1o+q1o');
        p1.Q2   = 0.5*(q2o+q2o');

        p1.name = sprintf('P%i:n=%i',p,n);
        p1.n    = n;
        p1.r    = r;

        plist{p,1} = p1;

    end

    %np = length(plist);

    %selectp     = 1:length(plist);
    nsel        = length(selectp);
    selplist    = plist(selectp);
    selns       = ns(selectp);
    selrs       = rs(selectp);

    %% Loop and output store
    for j = 1:nsel

        selProb     = selplist{j};
        prob.name   = selProb.name;

        clc;
        fprintf('Running Problem: %s \n',prob.name);

        g            = selProb.g;
        Q1           = selProb.Q1;
        Q2           = selProb.Q2;
        
        n            = size(g,1);
        prob.x0      = zeros(n,1);

        prob.n       = n;
        par.spIn     = speye(n);    
        prob.Pmat    = par.spIn;
        fncSB        = user_func_QuadStru(g,Q1,Q2);
        
        quadObj      = @(x)(x'*g+0.5*((x'*(Q1+Q2))*x));
        quadGrad     = @(x)(g+(Q1+Q2)*x);
    
        fncIP.objective = quadObj;
        fncIP.gradient  = quadGrad;

        fncLB        = {quadObj, quadGrad};

        %% Loop over Structured Solvers
        for s = 1:(nsolm2/2)
            
            scaling = mod(s,4);
            if scaling == 0; scaling = 4; end
            
            % Minus solver
            par.Scaling     = scaling;            
            outData{ir,j,s} = EX_CSBMSV3_F_MF_SCALE(prob,fncSB,par);
            
            % Plus solver
            parP.Scaling    = scaling;
            outData{ir,j,s+(nsolm2/2)} = EX_CSBP_SCALE_F(prob,fncSB,parP);
            
        end
        
        % External Solvers (IPOPT, LBFGSB)
        outData{ir,j,nsol-1} = EX_IPOPT_F(prob,fncIP,optsIP);
    
        outData{ir,j,nsol}   = EX_LBFGSB_F(prob,fncLB,optsLB);
        
        for k = 1:nsol

            % Averaging
            outIts(j,k)   = outIts(j,k)     + (outData{ir,j,k}.nIter/nruns);
            outObjs(j,k)  = outObjs(j,k)    + (outData{ir,j,k}.obj/nruns);
            outNgs(j,k)   = outNgs(j,k)     + (outData{ir,j,k}.ng/nruns);
            outTimes(j,k) = outTimes(j,k)   + (outData{ir,j,k}.ctime/nruns);

        end
        
        % Plot/store initialization figure
        name            = [PHI,'_','SIM',num2str(ir),'_',...
            'n',num2str(j)];
        
        dname           = [datapath,name];
        
        save(dname,'selProb');
        
    end
    
end

if saveFiles == 1
    
    dname           = [datapath,'experiments','_IV_',PHI];
    save(dname,'outData','outIts','outObjs','outNgs','outTimes',...
        'selplist','selns','selrs');

end

