% Experiment IVb:
% 06/28/17
% Comparison of 
% EIG-INF-G1-F-C,
% EIG-INF-G1
% L-BFGS-B (Stephen Becker).

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

paramsC = params; paramsC.c =4;


% Options for LBFGS-B
% Default setting for convergence.
opts = struct;
opts.m = 5;
opts.maxIts = 10000;
opts.maxTotalIts = 10000000;
opts.printEvery = Inf;
opts.factr=0;
 

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list.txt','r');  % cutest_list read file that contains CUTEr test problems

% Initialize storage for output information
numRuns=3; % 3
numAlgorithms=3;
numProblems=62;
ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);

p=1;
tline = fgets(fid);

while ischar(tline)     
    tline = fgets(fid);       
    
    if ~strcmp(tline(1),'%')  && ischar(tline)   
        
        eval(['!runcutest -p matlab -D ' tline]);
        %eval(['!runcuter --package mx --decode ' tline]);
        prob = cutest_setup();
        x0 = prob.x;
        
        n = size(x0,1);
        
        Ons = ones(n,1);
        l = -Inf.*Ons;
        u = +Inf.*Ons;
        
        opts.x0 = x0;
        opts.l = l;
        opts.u = u;
        
	% LBFGSB       
	s=1;        
	[ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LBFGS_B,x0,opts,numRuns);
        
        % EIG-INF-G1-F-C
        s=s+1;         
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G10_F_TERM_C,x0,paramsC,numRuns);

	% EIG-INF-G1
        s=s+1;       
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_TERM,x0,params,numRuns);
                        
        
        % Average CPU time
        if p==1 && numRuns > 2
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        else
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        end             
        
        cutest_terminate();
        
        p=p+1;            
    end
    
end

save('ResultsEXPERIMENT_VI','ex','numit','t_aver','numf','numg','params','tract');
