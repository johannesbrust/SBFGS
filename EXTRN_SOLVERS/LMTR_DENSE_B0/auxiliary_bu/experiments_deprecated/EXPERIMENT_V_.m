% Experiment V:
% Comparison of the values of gamma_perp as computed by EIG-INF-G1 and
% EIG-INF-G3. The ratio G1/G3 is computed.

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
sfil='TestResults';

% Initialize storage for output information
numRuns=10; % 3
numAlgorithms=3;
numProblems=62;

%{

ex=zeros(numProblems,numAlgorithms);
numf=zeros(numProblems,numAlgorithms);
numg=zeros(numProblems,numAlgorithms);
numit=zeros(numProblems,numAlgorithms);
tcpu=zeros(numProblems,numRuns,numAlgorithms);
t_aver=zeros(numProblems,numAlgorithms);
tract=zeros(numProblems,numAlgorithms);
numrst=zeros(numProblems,numAlgorithms);

%}

rats    = {};
rats_it = {};

probs   = [1 11 20 30 40 50];
np      = length(probs);
pnames  = {};

p       =1;
pc      = 1;
tline   = fgets(fid);

while ischar(tline)     
    tline = fgets(fid);       
    
    if ~strcmp(tline(1),'%')  && ischar(tline)   
        
        if probs(pc) == p
            
            pnames{pc} = tline;
            
        
            eval(['!runcutest -p matlab -D ' tline]);
            prob    = cutest_setup();
            x0      = prob.x;
            
            [ex,numf,numg,numit,tcpu,tract,numrst,numskip,rats{pc}, rats_it{pc}]=...
            runAlgorithm(@LMTR_EIG_inf_2_RAT,x0,params,numRuns);
                    
            cutest_terminate();
            
            if pc < np; pc = pc + 1; else pc = np; end;
        
        end
                
        p=p+1;            
    end
    
end

save('ResultsEXPERIMENT_V','ex','numit','numf','numg','params','rats', 'rats_it');

ylab    = '$ \gamma^{\perp}_k/ \hat{\gamma}^{\perp}_k  $';
xlab    = '$ k $';

for i = 1:np
    
    tit = pnames{i};
    plot(rats_it{i},rats{i}); 
            
    xlabel(xlab,'interpreter','latex');
    ylabel(ylab,'interpreter','latex');
    title(tit,'interpreter','latex');   
    
    name = [tit '.eps'];
    
    print('-dpsc2',tit);   %,-r200 
    
end


