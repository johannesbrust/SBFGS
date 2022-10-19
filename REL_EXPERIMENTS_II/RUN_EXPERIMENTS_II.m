%------------------------------ RUN_EXPERIMENTS_II -----------------------%
% RUN_EXPERIMENTS_II: Script to run the second experiments for the
% manuscript.
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
% Initial contributors: J.J.Brust, C.G.Petra, S.Leyffer.
%
%-------------------------------------------------------------------------%
%
% Initial version: J.B., 11/20/19, 
% 11/21/19, J.B., Very long test times. Preparing to patch data so 
% full memory methods only run onces.

clc;
clear;
close all;

addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));


tstart = tic;

probPath = 'fullProblemTable.txt';

%% Parameters
% Output containers with outputs for selected problems. 

dataPath        = fullfile('./data/');
figPath         = fullfile('./figs/');
saveFiles       = 1;

convTol         = 1e-6;
taumax          = 0;
legLocation     = 'SouthEast';

leg={           'L-S-BFGS-M1',...
                'L-S-BFGS-M2',...       
                'L-S-BFGS-M3',...
                'L-S-BFGS-M4',...
                'S-BFGS-M',...
                };
                
legP={          'L-S-BFGS-P1',...
                'L-S-BFGS-P2',...       
                'L-S-BFGS-P3',...
                'L-S-BFGS-P4',...
                'S-BFGS-P',...
                };
                
types.colors    = ['b' 'g' 'm' 'r' 'k'];

types.lines     = {'-.', '-', '-.','-','-.'};
                   
types.markers   = ['o' 'o' 'o' 'o' 'o'];


%% Experiments               
% Minus

% 11/21/19, J.B.
% Requires patching data. 

% m = 8;
 selectp = 1:61;
% EXPERIMENTS_II_A_FUNC(saveFiles,m,probPath,dataPath,nsol,selectp);
% %EXPERIMENTS_I_FUNC(saveFiles,PHI,figpath,datapath,nsol,nruns,ns,rscale);
% 
 m = 50;
 nsol = 4;
 EXPERIMENTS_II_A_FUNC(saveFiles,m,probPath,dataPath,nsol,selectp);

% patchData
fromData    = 'data/experiments_II_A_m8';
toData      = 'data/experiments_II_A_m50';
patchData(fromData,toData);

% Plus
% Load indices for small problems
% load('lsidx','sidx');
% %sidx = 37;
% m       = 8;
% nsol    = 5;
% EXPERIMENTS_II_B_FUNC(saveFiles,m,probPath,dataPath,nsol,sidx);
% 
% %sidx = 37;
% m       = 50;
% nsol    = 4;
% EXPERIMENTS_II_B_FUNC(saveFiles,m,probPath,dataPath,nsol,sidx);
% 
% % patchData
% fromData    = 'data/experiments_II_B_m8';
% toData      = 'data/experiments_II_B_m50';
% patchData(fromData,toData);

%% Plots
% Minus
selInd      = 1:5; % 1:4
%figPath     = fullfile('./figs/');

% 11/21/19, J.B.
% Requires patching data. 

m = 8;
expName     = ['II_A_m',num2str(m)];
GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);

m = 50;
expName     = ['II_A_m',num2str(m),'_patch'];
GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);
% 
% % Plus
% m = 8;
% expName     = ['II_B_m',num2str(m)];
% GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
%     legLocation,legP,types);
% 
% m = 50;
% expName     = ['II_B_m',num2str(m),'_patch'];
% GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
%     legLocation,legP,types);

tend = toc(tstart);