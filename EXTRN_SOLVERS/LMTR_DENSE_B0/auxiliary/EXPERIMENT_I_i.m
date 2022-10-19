%{
    Experiment I i:
    07/01/17, J.B

    Using unconstrained minimizer from multiple of identity matrix with
    two choices for gamma_perp:

    Choice 1: gamma_perp = c max_i gamma_i, c = 1,2,4,
    Choice 2: gamma_perp = lambda max_i gamma_i + (1-lambda) gamma, lambda
    = 1/4,2/4,3/4.

%}

% Experiment I:
% 06/27/17, J.B
% Comparison of large values of c for EIG-INF-G1 and EIG-INF-G1-F
% 
%           gamma_perp = c max_i gamma_i,
%           c = 1,2,4


clc;
clear;

parentpath = [cd(cd('..')) '/main'];
addpath(parentpath)

% Set input parameters:
params=struct;
params.m = 5;  % number of L-BFGS updates
params.gtol = 1e-5;  % exit if ||g||_2<gtol*max(1,||x||_2)
params.ranktol = 1e-7;  % tolerance for establishing rank of V
params.dflag = 0;  % display parameter, 1 if to display information per iteration
params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
params.ftol=1e-11;  % tolerance on relative function reduction

params.maxit = 10000;

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list_ex.txt','r');  % read file that contains CUTEr test problems

% Initialize storage for output information
numRuns=10; % 3
numAlgorithms=6;
numProblems=65;
ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);


numc =3;
parsC = cell(2*numc);


for i = 1:numc
   c = 2^((i-1));
   params.c = c;
   parsC{i} = params;
   
   params.c =1;
   params.lamc = i/4;
   parsC{i+numc} = params;   
end

%cl = parcluster('local');
%cl.parpool(4);

p=1;
tline = fgets(fid);
while ischar(tline)
    tline = fgets(fid);
    
    if ~strcmp(tline(1),'%')  && ischar(tline)
        
        eval(['!runcutest -p matlab -D ' tline]);
        %eval(['!runcuter --package mx --decode ' tline]);
        prob = cutest_setup();
        x0 = prob.x;
        
        % Multiple of identity for unconstrained minimizer
        % First three functions for c=1,2,4
        s=1; 
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % c=2
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % c=4
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % Last three functions for lambda=1/4,2/4,3/4
        % lambda=1/4
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C_L,x0,parsC{s},numRuns);
        
        % lambda=2/4
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C_L,x0,parsC{s},numRuns);
        
        % lambda=3/4
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C_L,x0,parsC{s},numRuns);
        
       
        % Average CPU time
        %for p = 1:numProblems
        if p==1 && numRuns > 2
            for si=1:numAlgorithms
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        else
            for si=1:numAlgorithms
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        end
        %end
        
        cutest_terminate();
        
        p=p+1;
        
    end
    
end
save('ResultsEXPERIMENT_I_i_10','ex','numit','t_aver','numf','numg','params');
