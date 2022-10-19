% Script to replot figures for paper using LaTeX Symbols.
% This script requires that the 'perf' function uses the argument pair
% 'Interpreter','Latex'.

clc;clear;

%% (1). Graph for comparison of dense initial matrix and multiple of identity:
load('ResultsEXPERIMENT_I');
leg={'$\left( \mathbf{P}, \infty \right) $',...
       '$\left( \mathbf{P}, \infty \right):\gamma^{\perp}_{k}$',...
       '$\left( \mathbf{P}, \infty \right):\hat{\gamma}^{\perp}_{k}$'};
         
indAlg = [1 2 3];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iter.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 time.eps;

close ALL;

%% (2). Graph for comparison of factorization methods.
load('ResultsEXPERIMENT_II');
leg={'$\left( \mathbf{P}, \infty \right)$:QR',...
       '$\left( \mathbf{P}, \infty \right)$:SVDI',...
       '$\left( \mathbf{P}, \infty \right)$:SVDII'};
   
indAlg = [1 2 3];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1);
print -dpsc2 iterSVD.eps;
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1);
print -dpsc2 timeSVD.eps;

close ALL;

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

close ALL;
 