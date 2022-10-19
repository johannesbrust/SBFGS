% Experiment IIIa:
% Comparison of small values of c for EIG-INF-G1
% 
%           gamma_perp = c max_i gamma_i,
%           c = 1,1/2,1/4,1/8


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
fid =fopen('cutest_list.txt','r');  % read file that contains CUTEr test problems

% Initialize storage for output information
numRuns=10; % 3
numAlgorithms=4;
numProblems=62;
ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);


numc =4;
parsC = cell(numc);

for i = 1:numc
   c = 2^(-(i-1));
   params.c = c;
   parsC{i} = params;
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
        
        % c=1
        s=1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % c=1/2
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % c=1/4
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        % c=1/8
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_C,x0,parsC{s},numRuns);
        
        % Average CPU time
        if p==1 && numRuns > 2
            for si=1:numc
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        else
            for si=1:numc
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        end
        
        cutest_terminate();
        
        p=p+1;
        
    end
    
end

save('ResultsEXPERIMENT_IIIa','ex','numit','t_aver','numf','numg','params');