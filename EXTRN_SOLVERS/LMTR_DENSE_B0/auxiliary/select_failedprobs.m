% Script to identify the problems on which all solvers, in all experiments
% failed.
% J.B., 09/22/17

% Experiment 1 b. (SCALING)
load('ResultsEXPERIMENT_I_i_10'); indAlg = [1 2 5 4];
ex = ex(:,indAlg);
f1b = sum(ex > 0, 2) == 0;

% Experiment 1 a.
load('ResultsEXPERIMENT_I_ii_10'); indAlg = [1 2 8 7];
ex = ex(:,indAlg);
f1a = sum(ex > 0, 2) == 0;

% Experiment 2. (DENSE)
load('ResultsEXPERIMENT_II_i'); indAlg = [1 2];
ex = ex(:,indAlg);
f2 = sum(ex > 0, 2) == 0;

% Experiment 3. (SVD)
load('ResultsEXPERIMENT_III_i'); indAlg = [1 2 3];
ex = ex(:,indAlg);
f3 = sum(ex > 0, 2) == 0;

% Experiment 4. (LBFGSB)
load('ResultsEXPERIMENT_VI_ii_10'); indAlg = [2 1];
ex = ex(:,indAlg);
f4 = sum(ex > 0, 2) == 0;

% Experiment 5. (TR)
load('ResultsEXPERIMENT_V_ii_10'); indAlg = [2 1];
ex = ex(:,indAlg);
f5 = sum(ex > 0, 2) == 0;


% Experiment 6. (ORIGINAL)
load('ResultsEXPERIMENT_IV_ii_10'); indAlg = [1 2];
ex = ex(:,indAlg);
f6 = sum(ex > 0, 2) == 0;

F = [f1a f1b f2 f3 f4 f5 f6]


