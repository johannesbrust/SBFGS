%---------------------- Test patchData -----------------------------------%
%
% Script to test the patchData function
%
%-------------------------------------------------------------------------%
% 11/21/19, J.B., Initial Version

fromData    = 'data/experiments_II_B_m8';
toData      = 'data/experiments_II_B_m50';

patchData(fromData,toData);

% Expect data 'data/experiments_II_B_m50_patch'; 

lpatch = load('data/experiments_II_B_m50_patch');