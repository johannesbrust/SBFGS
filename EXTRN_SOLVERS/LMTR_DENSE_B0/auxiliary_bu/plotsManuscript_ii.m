% Script to replot figures for paper using LaTeX Symbols.
% This script requires that the 'perf' function uses the argument pair
% 'Interpreter','Latex'.

% 06/27/17, Latest setup of experiments.
% 06/28/17 Plots according to description in "Comments (part 3) on version
% of 27.05.2017 " (from DenseB0_SIOPT/notes)
% 07/05/17, Updates to comparison plots.
% 07/09/17 Updates to allow for pairwise comparisons.
% 07/12/17 Replotting results in reference to email 7/10/17
% 07/21/17 Replotting results for comparisons
% 09/19/17 Selection of plots for Experiment 4.
% 09/20/17 Legend to Experiment 5., for consistency with manuscript 
% 09/22/17 Replotting to exclude graphs where are all solvers failed
clc;clear;

whichcomp =3;   % Flag to determine comparison
                % := 1, compare three solvers,
                % := 2, compare other \hat{p}_u,
                % := 3, compare other p_u.

%% (1 i). Comparison for p_u using c = 1,2,4, lambda = 1/4,2/4,3/4
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_i_10'); % 'ResultsEXPERIMENT_I_i'
leg={'$\hat{\mathbf{B}}_0(1,1)$',...
       '$\hat{\mathbf{B}}_0(2,1)$',...
       '$\hat{\mathbf{B}}_0(4,1)$',...
       '$\hat{\mathbf{B}}_0(1, \frac{1}{4})$',...
       '$\hat{\mathbf{B}}_0(1, \frac{1}{2})$',...
       '$\hat{\mathbf{B}}_0(1, \frac{3}{4})$'};
            
%indAlg = [1 2 3 4 5 6]; % 4 5 6
indAlg  = [1 2 5 4];
[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);

% Set line types
types.colors = ['b' 'g' 'r' 'k'];
types.lines = {'-' '--' ':' '-.'};
types.markers = ['s' 'o'];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex1_b.eps; savefig('iter_ex1_b');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex1_b.eps; savefig('time_ex1_b');

close ALL;

%% (1 ii). Comparison for \hat{p}_u using c = 1,2,4,8,10,100 lambda = 1/4,2/4,3/4
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_ii_10'); % ResultsEXPERIMENT_I_ii
leg={  '$\hat{\mathbf{B}}_0(1,1)^*$',...
       '$\hat{\mathbf{B}}_0(2,1)^*$',...
       '$\hat{\mathbf{B}}_0(4,1)^*$',...
       '$\hat{\mathbf{B}}_0(8,1)^*$',...
       '$\hat{\mathbf{B}}_0(10,1)^*$',...
       '$\hat{\mathbf{B}}_0(100,1)^*$',...
       '$\hat{\mathbf{B}}_0(1, \frac{1}{4})^*$',...
       '$\hat{\mathbf{B}}_0(1, \frac{1}{2})^*$',...
       '$\hat{\mathbf{B}}_0(1, \frac{3}{4})^*$'};
            
%indAlg = [1 2 3 4 5 6 7 8 9]; % 4 5 6

indAlg = [1 2 8 7];

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);

% Set line types
types.colors = ['b' 'g' 'r' 'k'];
types.lines = {'-' '--' ':' '-.'};
types.markers = ['s' 'o'];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex1_a.eps; savefig('iter_ex1_a');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex1_a.eps; savefig('time_ex1_a');
close ALL;

%% (1 iii). Comparison of the best of (1 i) and (1 ii)
% Data also available for same experiments without the formula
load('ResultsEXPERIMENT_I_i_10');
indi = [1 5];
exi = ex(:,indi); numiti = numit(:,indi); t_averi = t_aver(:,indi);

load('ResultsEXPERIMENT_I_ii_10');
indii = [1 8];
exii = ex(:,indii); numitii = numit(:,indii); t_averii = t_aver(:,indii);

leg={  '$\hat{\mathbf{B}}_0(1,1)$',...
       '$\hat{\mathbf{B}}_0(1,\frac{1}{2})$',...
       '$\hat{\mathbf{B}}_0(1,1)^*$',...
       '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$'};
            
%indAlg = [1 2 3 4]; % 4 5 6
%indAlg = [2 4];

indAlg = [1 4];

types.colors = ['b' 'r'];
types.lines = {'-' ':'};
types.markers = ['s' 'o'];


ex = [exi exii]; numit = [numiti numitii]; t_aver =[t_averi t_averii];

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1, types);
print -dpsc2 iter_ex2.eps; savefig('iter_ex2');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1, types);
print -dpsc2 time_ex2.eps; savefig('time_ex2');

close ALL;

%% Experiment II-Comparison of best solver from experiments a,b
load('ResultsEXPERIMENT_II_i');

leg={  '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$',...
       '$\hat{\mathbf{B}}_0(1,\frac{1}{2})$'};
            
%indAlg = [1 2 3 4]; % 4 5 6
%indAlg = [2 4];

indAlg = [1 2];

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);

types.colors = ['r' 'b'];
types.lines = {':' '-'};
types.markers = ['s' 'o'];


%ex = [exi exii]; numit = [numiti numitii]; t_aver =[t_averi t_averii];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1, types);
print -dpsc2 iter_ex2_i.eps; savefig('iter_ex2_i');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1, types);
print -dpsc2 time_ex2_i.eps; savefig('time_ex2_i');

close ALL;


%% Experiment III-Comparison of best solver from experiment II with SVD Solvers
load('ResultsEXPERIMENT_III_i');
% load('ResultsEXPERIMENT_III_i_a');

leg={  '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$-QR',...
       '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$-SVD I',...
       '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$-SVD II'};
            
%indAlg = [1 2 3 4]; % 4 5 6
%indAlg = [2 4];

indAlg = [1 2 3];

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);

types.colors = ['r' 'b' 'g'];
types.lines = {':' '-.' '-'};
types.markers = ['s' 'o'];


%ex = [exi exii]; numit = [numiti numitii]; t_aver =[t_averi t_averii];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1, types);
print -dpsc2 iter_ex3_i.eps; savefig('iter_ex3_i');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1, types);
print -dpsc2 time_ex3_i.eps; savefig('time_ex3_i');

close ALL;



%% (4 i). Comparison of best of (1 i) and (1 ii) to original LMTR method
%load('ResultsEXPERIMENT_IV_i_10');
load('ResultsEXPERIMENT_IV_ii_10');

leg={  '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$',...
       '$\mathbf{B}_0(\gamma_k)$'};

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);
   
types.colors = ['r' 'k'];
types.lines = {':' '-.'};
types.markers = ['s' 'o'];

indAlg = [1 2];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex4_i_10.eps; savefig('iter_ex4_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex4_i_10.eps; savefig('time_ex4_i_10');

% Plot results for the iterations where steplenght=1 is rejected often. In particular, 
% where the step-size one was 
% rejected in, at least, 30% of iterations

indProb=find(tract(:,1)./numit(:,1)>=0.30);
perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex4_i_10_sel.eps; savefig('iter_ex4_i_10_sel');
perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex4_i_10_sel.eps; savefig('time_ex4_i_10_sel');

close ALL;
   

%% (5 i). Comparison of best of (1 i) and (1 ii) to L-BFGS-B
%load('ResultsEXPERIMENT_V_i_10');
load('ResultsEXPERIMENT_V_ii_10');

 leg={ 'L-BFGS-B',... 
     '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$'};
         
[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);
 
types.colors = ['k' 'r'];
types.lines = {'-' ':'};
types.markers = ['s' 'o'];

indAlg = [2 1];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex5_i_10.eps; savefig('iter_ex5_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex5_i_10.eps; savefig('time_ex5_i_10');

close ALL;

%% (6 i). Comparison of best of (1 i) and (1 ii) to L-BFGS-TR
%load('ResultsEXPERIMENT_VI_i_10');
load('ResultsEXPERIMENT_VI_ii_10');

 leg={ 'L-BFGS-TR',... 
     '$\hat{\mathbf{B}}_0(1,\frac{1}{2})^*$'};   

[n,m]   = dim(ex); idxr = [1:7;11:n];
ex      = ex(idxr,:);
numit   = numit(idxr,:);
t_aver  = t_aver(idxr,:);
 
types.colors = ['k' 'r'];
types.lines = {'-' ':'};
types.markers = ['s' 'o'];

indAlg = [2 1];

perf(ex(:,indAlg),numit(:,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex6_i_10.eps; savefig('iter_ex6_i_10');
perf(ex(:,indAlg),t_aver(:,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex6_i_10.eps; savefig('time_ex6_i_10');


% Plot results for the iterations where steplenght=1 is rejected often. In particular, 
% where the step-size one was 
% rejected in, at least, 30% of iterations

indProb=find(tract(:,1)./numit(:,1)>=0.30);
perf(ex(indProb,indAlg),numit(indProb,indAlg),leg(indAlg),1,types);
print -dpsc2 iter_ex6_i_10_sel.eps; savefig('iter_ex6_i_10_sel');
perf(ex(indProb,indAlg),t_aver(indProb,indAlg),leg(indAlg),1,types);
print -dpsc2 time_ex6_i_10_sel.eps; savefig('time_ex6_i_10_sel');


close ALL;

 