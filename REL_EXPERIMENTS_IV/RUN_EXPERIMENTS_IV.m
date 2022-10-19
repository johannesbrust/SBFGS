%------------------------------ RUN_EXPERIMENTS_IV -----------------------%
% RUN_EXPERIMENTS_IV: Script to run the fourth experiment of the
% manuscript. This is about comparing to external solvers.
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
% Initial version: J.B., 11/20/19, Test small ~7min (4.146856330870000e+02 secs.)
% 12/03/19, J.B., Test large clustered eigenvalues

clc;
clear;
close all;

addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));
addpath(genpath('../EXTRN_SOLVERS/IPOPT'));
addpath(genpath('../EXTRN_SOLVERS/LBFGSB'));

tstart = tic;

%% Parameters
% Output containers with outputs for selected problems. 
nsol            = 10; %5
nruns           = 2;

ns              = [100;200;300;400;500;600;700]; % [100;200]% ;300;400
rscale          = 10;

datapath        = fullfile('./data/');
figpath         = fullfile('./figs/');
saveFiles       = 1;

convTol         = 5*1e-6;
taumax          = 0;
legLocation     = 'SouthEast';

leg={           'L-S-BFGS-M1',...
                'L-S-BFGS-M2',...       
                'L-S-BFGS-M3',...
                'L-S-BFGS-M4',...
                'L-S-BFGS-P1',...
                'L-S-BFGS-P2',...       
                'L-S-BFGS-P3',...
                'L-S-BFGS-P4',...
                'IPOPT',...
                'L-BFGS-B'
                };
                
types.colors    = ['b' 'g' 'm' 'r'...
                   'b' 'g' 'm' 'r' 'k' 'y']; %'k' 'y'
types.lines     = {'-.', '-', '-.','-'...
                   '-.', '-', '-.','-','-.','-'}; %'-',   '-'
types.markers   = ['o' 'o' 'o' 'o'...
                   'o' 'o' 'o' 'o' 'o' 'o']; %'s' 's'

%% Experiments               
selectp = 1:length(ns);
PHI     = 'SM'; % SM
EXPERIMENTS_IV_FUNC(saveFiles,PHI,datapath,nsol,nruns,ns,rscale,selectp);

%PHI = 'LRG';
%EXPERIMENTS_IV_FUNC(saveFiles,PHI,datapath,nsol,nruns,ns,rscale,selectp);

%% Minus Plots
% First small, then large eigenvalue clustering
selInd      = [1:4,9,10]; % 1:4
figPath     = fullfile('./figs/Minus');

expName     = 'IV_SM'; % IV_SM
GEN_PLOTS_FUNC(expName,datapath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);

% expName     = 'I_LRG'; 
% GEN_PLOTS_FUNC(expName,datapath,figPath,selInd,convTol,taumax,...
%     legLocation,leg,types);

%% Plus Plots
selInd      = [5:8,9,10];
figPath     = fullfile('./figs/Plus'); 

expName     ='IV_SM'; % IV_SM
GEN_PLOTS_FUNC(expName,datapath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);

% expName     = 'I_LRG'; 
% GEN_PLOTS_FUNC(expName,datapath,figPath,selInd,convTol,taumax,...
%     legLocation,leg,types);

tend = toc(tstart);