%------------------------------ RUN_EXPERIMENTS_III ----------------------%
% RUN_EXPERIMENTS_III: Script to run the application experiments for the
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
% Initial version: J.B., 11/20/19, 

clc;
clear;
close all;

addpath(genpath('../ALG_COMPACT'));
addpath(genpath('../INTERFACE'));
addpath(genpath('../MISC'));


tstart = tic;

LIBSVMPath = '/Users/johannesbrust/Dropbox/codes/LIBSVM/';

%% Parameters
% Output containers with outputs for selected problems. 
 %5

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
                'BFGS',...
                };
                
% legP={          'L-S-BFGS-P1',...
%                 'L-S-BFGS-P2',...       
%                 'L-S-BFGS-P3',...
%                 'L-S-BFGS-P4',...
%                 'S-BFGS-P',...
%                 };
                
types.colors    = ['b' 'g' 'm' 'r' 'k'];

types.lines     = {'-.', '-', '-.','-','-.'};
                   
types.markers   = ['o' 'o' 'o' 'o' 'o'];


%% Experiments               
% LIBSVM Minus
nsol = 4;
selectp = 1:10;
EXPERIMENTS_III_A_FUNC(saveFiles,dataPath,LIBSVMPath,nsol,selectp);

% PDE Minus and BFGS
N   = [20,30,40,50,60,70,80,90,100];
xl  = 0;
xu  = 1;
yl  = 0;
yu  = 1;
a   = 0.0;
b   = 0.0;
sig = 1.0; % 1.0, 5.0

nsol = 5;
EXPERIMENTS_III_B_FUNC(saveFiles,dataPath,N,xl,xu,yl,yu,a,b,sig,nsol);

%% Plots
% LIBSVM

selInd = 1:4;
expName     = 'III_A';
GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);

% PDE

selInd = 1:5;
expName     = 'III_B';
GEN_PLOTS_FUNC(expName,dataPath,figPath,selInd,convTol,taumax,...
    legLocation,leg,types);

tend = toc(tstart);