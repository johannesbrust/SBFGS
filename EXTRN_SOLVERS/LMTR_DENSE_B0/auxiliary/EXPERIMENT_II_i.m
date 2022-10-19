% Experiment II:
% 06/27/17, J.B
% Comparison of using a formula for computing the unconstrained minimizer
% EIG-INF-G1 (Proposed method),
% EIG-INF-G10-F (Dense method, max gamma_k).

% Updated: 07/13/17 Newly proposed experiment II, i.e
%           (hat{B}_0, lambda=1/2)* and hat{B}_0

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

paramsC = params;
paramsC.lamc = 1/2;
paramsC.c = 1;



params.maxit = 10000;

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list_ex.txt','r');  % cutest_list, read file that contains CUTEr test problems
%sfil='TestResults';

% Initialize storage for output information
numRuns=10; % 3
numAlgorithms=2;
numProblems=65;
ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);
numidec=zeros(numProblems,numAlgorithms);

p=1;
tline = fgets(fid);
while ischar(tline)     
    tline = fgets(fid);       
    
    if ~strcmp(tline(1),'%')  && ischar(tline)   
        
        eval(['!runcutest -p matlab -D ' tline]);
        %eval(['!runcuter --package mx --decode ' tline]);
        prob = cutest_setup();
        x0 = prob.x; 
                
        % Formula
         s=1;        
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s),numidec(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G10_F_C_L,x0,paramsC,numRuns);
        
        % Multiple of identity
        s=s+1;        
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s),numidec(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1,x0,params,numRuns);
        
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
    %save(sfil);
end

save('ResultsEXPERIMENT_II_i','ex','numit','t_aver','numf','numg','params');
