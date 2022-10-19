% This is a script to minimize functions from the CUTEst test set using 
% the Eig(inf,2) and and L-BFGS method. The intend is to reproduce Fig.2 of
% the paper below.
% J.B 10/13, This script is modified for the CUTEst test suite.
%
% For details, see:
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Link�ping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359

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

% Set maximum iterations, J.B 20/10
params.maxit = 1000;

% J.B 02/11/16
% Additional options for lbfgsb

opts.m          = params.m;
opts.factr      = 0;
opts.pgtol      =1e-5;
opts.maxIts     = 1000;
opts.printEvery = 0;
params.opts     = opts;

CUTEst_init  % initialize CUTEr, see appropriate documentation 
fid =fopen('cutest_list.txt','r');  % read file that contains CUTEr test problems
sfil='TestResults';

% 02/11/16 Adding the L-BFGS-B algorithm by Nocedal et. al
leg={'L-BFGS-MTBT','L-BFGS','L-BFGS-B', 'EIG(\infty,2)','m=10','m=15','EIG-MS','EIG-MS(2,2)','BWX-MS','D-DOGL'};

params10=params;
params10.m=10;
params15=params;
params15.m=15;

% Initialize storage for output information
numRuns=1;
numAlgorithms=length(leg);
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
        
        onesv           = ones(size(x0,1),1);
        params.l        = -inf.*onesv;
        params.u        = inf.*onesv;
        params.opts.x0  = x0;
                
         s=1;        
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LBFGS_MTBT,x0,params,numRuns);
        
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LBFGS_TR,x0,params,numRuns);
        
        s=s+1;
        [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
            runAlgorithm(@LBFGS_B,x0,params,numRuns);
        
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2,x0,params,numRuns);  
        
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2,x0,params10,numRuns); 
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_inf_2,x0,params15,numRuns);
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_MS,x0,params,numRuns);        
%                 
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s),numrst(p,s)]=...
%             runAlgorithm(@LMTR_EIG_MS_2_2,x0,params,numRuns);        
%         
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_BWX_MS,x0,params,numRuns);                                        
%                 
%         s=s+1;
%         [ex(p,s),numf(p,s),numg(p,s),numit(p,s),tcpu(s,:,p),tract(p,s)]=...
%             runAlgorithm(@LMTR_DDOGL,x0,params,numRuns);                
        
        % Average CPU time
        if p==1
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,3:numRuns,p))/(numRuns-2);
            end
        else
            for si=1:s
                t_aver(p,si) = sum(tcpu(si,2:numRuns,p))/(numRuns-1);
            end
        end             
        
        p=p+1;            
    end
    save(sfil);
end


%% Plot performance profiles

% Fig.1: EIG(inf,2) for m=5,10,15
% indAlg=[3,4,5];
% leg(3)={'m=5'};
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));

% Fig.2: EIG(inf,2) and L-BFGS
indAlg=[3,2];
%leg(3)={'EIG(\infty,2)'};
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));

% Fig.3: EIG(inf,2) and L-BFGS in those problems where the step-size one was 
% rejected in, at least, 30% of iterations
% indProb=find(tract(:,2)./numit(:,2)>=0.30);
% perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg));
% perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg));
% 
% % Fig.4: EIG(inf,2), EIG-MS and EIG-MS(2,2)
% indAlg=[3,6,7];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
% 
% % Fig.5: EIG(inf,2) and BWX-MS
% indAlg=[3,8];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
% 
% % Fig.6: EIG(inf,2) and D-DOGL
% indAlg=[3,9];
% perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg));
% perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg));
