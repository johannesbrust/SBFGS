% Experiment III:
% 06/28/17, J.B 
% 07/16/17
% QR 
%   - EIG-INF-G1-F-Lambda-1/2
% SVDI (SVD based),
% SVDII (SVD based, using the approach described by X.Lu 
%        "A Study of the Limited Memory SR1 Method in Practice").

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

paramsC = params; paramsC.c =1; paramsC.lamc=1/2; % Scaling to the dense parameter

paramsC.maxit = 10000; 

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list_ex.txt','r');  % read file that contains CUTEr test problems
sfil='TestResults';

% Initialize storage for output information
numRuns=10; % 3
numAlgorithms=3;
numProblems=65;
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
        prob = cutest_setup();
        x0 = prob.x; 
                
        % EIG-INF-G1-F-C
         s=1;        
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G10_F_C_L,x0,paramsC,numRuns);
        
        
        % SVD I
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_SVD1_C_L,x0,paramsC,numRuns);
                
        % SVDII
        s = s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LMTR_EIG_inf_2_DENSE_G1_SVD2_C_L,x0,paramsC,numRuns);
        
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

save('ResultsEXPERIMENT_III_i','ex','numit','t_aver','numf','numg','params');
