function [figout] = GEN_PLOTS_FUNC(expName, dataPath, figPath, selInd, convTol,...
    taumax, legLocation, leg, types)
% GEN_PLOTS_FUNC: Function to plot data from SBFGS numerical experiments.
%
% Extension to structured BFGS formulas from Petra, Chiang,
% Anitescu 2018. The extended methods use compact representations of the
% quasi-Newton formulas. The compact forms are
%
% B = B0 - Psi M Psi',
%
% where typically B0 = gamma.I (n x n) multiple of identity initial matrix,
% Psi (n x 2m), M (2m x 2m) are small low rank updates. 
%
% Initial contributors: J.J.Brust, C.G.Petra, S.Leyffer, M.Anitescu.
%
% INPUTS:
% expName   := Experiment ending for data loading
% dataPath  := Data path
% figPath   := Figure path
% selInd    := Selection of indices
% convTol   := Convergence tolerance
% taumax    := Possibly to truncate performence plots
% legLocation := Legend location
% leg       := Legend
% types     := Line specifiers
%
% OUTPUS:
% figout    := Last figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11/19/19, J.B

   
expType = 'ITER';
% Plot and store ITER
EXPERIMENTS_PLOT(expName, expType, figPath, dataPath,...
    convTol, leg, types, taumax, selInd, legLocation);
title('ITER');

% Plot and store TIME
expType = 'TIME';
figout = EXPERIMENTS_PLOT(expName, expType, figPath, dataPath,...
    convTol, leg, types, taumax, selInd, legLocation);
title('TIME');