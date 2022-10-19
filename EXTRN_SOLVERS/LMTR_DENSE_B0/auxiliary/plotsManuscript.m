% Script to replot figures for paper using LaTeX Symbols.
% This script requires that the 'perf' function uses the argument pair
% 'Interpreter','Latex'.

% 06/27/17, Latest setup of experiments.
% 06/28/17 Plots according to description in "Comments (part 3) on version
% of 27.05.2017 " (from DenseB0_SIOPT/notes)
% 07/05/17, Updates to comparison plots.
% 07/09/17 Updates to allow for pairwise comparisons.
clc;clear;

whichcomp =3;   % Flag to determine comparison
                % := 1, compare three solvers,
                % := 2, compare other \hat{p}_u,
                % := 3, compare other p_u.

%% (1). Comparison for \hat{p}_u using \gamma_perp = c max \gamma_i
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I');
leg={'$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,2\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,2\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)$'};
            
indAlg = [1 2 3]; % 4 5 6

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex1.eps; savefig('iter_ex1');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex1.eps; savefig('time_ex1');

% Plots for internal comparison

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex1_.eps; savefig('iter_ex1_');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex1_.eps; savefig('time_ex1_');

close ALL;

%% (1 i). Comparison for p_u using c = 1,2,4, lambda = 1/4,2/4,3/4
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_i_10'); % 'ResultsEXPERIMENT_I_i'
leg={'$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,2\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 1/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 3/4)$'};
            
indAlg = [1 2 3 4 5 6]; % 4 5 6

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex1_i.eps; savefig('iter_ex1_i');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex1_i.eps; savefig('time_ex1_i');

close ALL;

%% (1 ii). Comparison for \hat{p}_u using c = 1,2,4,8,10,100 lambda = 1/4,2/4,3/4
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_ii_10'); % ResultsEXPERIMENT_I_ii
leg={  '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,2\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,8\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,10\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,100\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 1/4)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 3/4)^*$'};
            
indAlg = [1 2 3 4 5 6 7 8 9]; % 4 5 6

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex1_ii.eps; savefig('iter_ex1_ii');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex1_ii.eps; savefig('time_ex1_ii');

close ALL;

%% (1 iii). Comparison of the best of (1 i) and (1 ii)
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_i_10');
indi = [1 5];
exi = ex(:,indi); numiti = numit(:,indi); t_averi = t_aver(:,indi);

load('ResultsEXPERIMENT_I_ii_10');
indii = [1 8];
exii = ex(:,indii); numitii = numit(:,indii); t_averii = t_aver(:,indii);

leg={  '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 3/4)^*$'};
            
indAlg = [1 2 3 4]; % 4 5 6
%indAlg = [2 4];

ex = [exi exii]; numit = [numiti numitii]; t_aver =[t_averi t_averii];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex1_iii.eps; savefig('iter_ex1_iii');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex1_iii.eps; savefig('time_ex1_iii');

close ALL;


%% (4 i). Comparison of best of (1 i) and (1 ii) to original LMTR method
load('ResultsEXPERIMENT_IV_i_10');

leg={  '$\mathbf{B}_0(\gamma_k)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

if whichcomp == 1
    indAlg = 1:3;
elseif whichcomp == 2
        indAlg = [1 3];   
else
    indAlg = [1 2];
end

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex4_i_10.eps; savefig('iter_ex4_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex4_i_10.eps; savefig('time_ex4_i_10');

close ALL;
   

%% (5 i). Comparison of best of (1 i) and (1 ii) to L-BFGS-B
load('ResultsEXPERIMENT_V_i_10');

leg={  'L-BFGS-B',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

if whichcomp == 1
    indAlg = 1:3;
elseif whichcomp == 2
        indAlg = [1 3];   
else
    indAlg = [1 2];
end
   
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex5_i_10.eps; savefig('iter_ex5_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex5_i_10.eps; savefig('time_ex5_i_10');

close ALL;

%% (6 i). Comparison of best of (1 i) and (1 ii) to L-BFGS-TR
load('ResultsEXPERIMENT_VI_i_10');

leg={  'L-BFGS-TR',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

if whichcomp == 1
    indAlg = 1:3;
elseif whichcomp == 2
        indAlg = [1 3];   
else
    indAlg = [1 2];
end
   
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex6_i_10.eps; savefig('iter_ex6_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex6_i_10.eps; savefig('time_ex6_i_10');


% Plot results for the iterations where steplenght=1 is rejected often. In particular, 
% where the step-size one was 
% rejected in, at least, 30% of iterations

indProb=find(tract(:,1)./numit(:,1)>=0.30);
perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex6_i_10_sel.eps; savefig('iter_ex6_i_10_sel');
perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg),1);
print -dpsc2 time_ex6_i_10_sel.eps; savefig('time_ex6_i_10_sel');


close ALL;

   
%{
%% (2). Comparison of \hat{p}_u to p_u.
% Data available in 'ResultsEXPERIMENT_I'
load('ResultsEXPERIMENT_I');
leg={'$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)^*$','','',...
       '$$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};
   
indAlg = [1 4];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex2.eps; savefig('iter_ex2.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex2.eps; savefig('time_ex2.fig');

close ALL;

%% (3). Comparison of SVD approaches with best choice of experiments (1) and (2).
load('ResultsEXPERIMENT_III');
leg={'$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$-SVDI',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$-SVDII'};

% Plots for internal comparison
indAlg = [1 2 3 4];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex3_.eps; savefig('iter_ex3_.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex3_.eps; savefig('time_ex3_.fig');

% Plots according to document, best choice 4 \gamma^max
indAlg = [1 3 4];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex3.eps; savefig('iter_ex3.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex3.eps; savefig('time_ex3.fig');

close ALL;

%% (4). Comparison to LMTR with best choice of experiments (1) and (2).
load('ResultsEXPERIMENT_IV');
leg={'$\mathbf{B}_0(\gamma_k)$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

% Plots for internal comparison
    
indAlg = [1 2 3];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex4_.eps; savefig('iter_ex4_.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex4_.eps; savefig('time_ex4_.fig');

% Plots according to document
indAlg = [1 2];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex4.eps; savefig('iter_ex4.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex4.eps; savefig('time_ex4.fig');

%% (5). Comparison of best choice to LBFGSB.
load('ResultsEXPERIMENT_VI'); % Order of experiments unaligned
leg={'L-BFGS-B',...
        '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

% Plots for internal comparison
    
indAlg = [1 2 3];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex5_.eps; savefig('iter_ex5_.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex5_.eps; savefig('time_ex5_.fig');

% Plots according to document
indAlg = [1 2];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex5.eps; savefig('iter_ex5.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex5.eps; savefig('time_ex5.fig');

%% (6). Comparison of best choice to LBFGS-TR.
load('ResultsEXPERIMENT_V'); % Order of experiments unaligned
leg={'L-BFGS-TR',...
        '$\hat{\mathbf{B}}_0(\gamma_k,4\gamma^{\perp}_k)^*$',...
        '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};

% Plots for internal comparison
    
indAlg = [1 2 3];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex6_.eps; savefig('iter_ex6_.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex6_.eps; savefig('time_ex6_.fig');

% Plots according to document
indAlg = [1 2];
perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex6.eps; savefig('iter_ex6.fig');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex6.eps; savefig('time_ex6.fig');

%% (7). Comparison for flexible parameter.
load('ResultsEXPERIMENT_VII');
leg={'$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 1/4)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 1/2)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 3/4)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)^*$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 1/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 2/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\hat{\gamma}^{\perp}_k, \lambda = 3/4)$',...
       '$\hat{\mathbf{B}}_0(\gamma_k,\gamma^{\perp}_k)$'};
         
indAlg = [1 2 3 4 5 6 7 8];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter_ex7.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time_ex7.eps;

close ALL;

%}

%{
%% (3a). Graph for comparison of small values for c.
load('ResultsEXPERIMENT_IIIa');
leg={ '$\left( \mathbf{P}, \infty \right):c_{-1}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{-2}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{-3}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{-4}\gamma^{\perp}_{k}$'};
       
indAlg = [1 2 3 4];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterCSmall.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeCSmall.eps;

close ALL;

%% (3b). Graph for comparison of small values for c.
load('ResultsEXPERIMENT_IIIb');
leg={ '$\left( \mathbf{P}, \infty \right):c_{1}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{2}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{3}\gamma^{\perp}_{k}$',...
      '$\left( \mathbf{P}, \infty \right):c_{4}\gamma^{\perp}_{k}$'};
       
indAlg = [1 2 3 4];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterCLarge.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeCLarge.eps;

close ALL;


%% (4a). Graph for comparison of dense trust-region method with Oleg's line-search method.
load('ResultsEXPERIMENT_IVa');
leg={ '$\left( \mathbf{P}, \infty \right):\gamma^{\perp}_{k}$',...
       'L-BFGS-TR'};
       
indAlg = [1 2];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterLBFGSTR.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeLBFGSTR.eps;

% Plot results for the iterations where steplenght=1 is rejected often. In particular, 
% where the step-size one was 
% rejected in, at least, 30% of iterations

indProb=find(tract(:,2)./numit(:,2)>=0.30);
perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg),1);
print -dpsc2 iterLBFGSTR_SEL.eps;
perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg),1);
print -dpsc2 timeLBFGSTR_SEL.eps;


close ALL;

%% (4b). Graph for comparison of dense trust-region method with L-BFGS-B line-search method.
load('ResultsEXPERIMENT_IVb');
leg={ '$\left( \mathbf{P}, \infty \right):\gamma^{\perp}_{k}$',...
       'L-BFGS-B'};
       
indAlg = [1 2];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterLBFGSB.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeLBFGSB.eps;

close ALL;

%% (2 New) Experiments with using a formula for the dense initial matrix.
load('ResultsEXPERIMENT_II_F');
leg={'$\left( \mathbf{P}, \infty \right):\gamma^{\perp}_{k}$:F',...
     '$\left( \mathbf{P}, \infty \right):\hat{\gamma}^{\perp}_{k}$:F',...
     '$\left( \mathbf{P}, \infty \right):\gamma_{k}$'};
 
 indAlg = [1 2 3];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterF.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeF.eps;

%}

close ALL;
 