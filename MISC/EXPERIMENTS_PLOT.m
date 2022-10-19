function out = EXPERIMENTS_PLOT(expName, expType, figPath, dataPath,...
    convTol, leg, types, taumax, indAlg, legLocation)
% EXPERIMENTS_PLOT: Plots for Compact-structured BFGS experiments.
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
% An initial report on the compact representations is in DOCS/report_compact_sQN_011819 
%
% Function to plot results. 
%
% INPUTS:
% expName   := Experiment name. Values are {'I_A','I_B',...}. Creates file name:
% EXPERIMENT_expName_outType.pdf
% expType   := Experiment type. Values are {'ITER','TIME'}
% figPath   := Figure store path
% dataPath  := Data load path
% convTol   := Classify problem as solved if norm(g) <= convTol
% leg       := Plot legend
% types     := Line types for plot
% taumax    := Parameter to allow for a max. tau value
% indAlg    := Which algorithms (index algorithms)
% legLocation:= Legend location
%
% OUTPUTS:
% out       := Plot handel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11/15/19, J.B., Initial version 
% 11/19/19, J.B., Modification for algorithm indices

% Plotting parameters 

%figpath = './figs';
figname     = ['EXPERIMENT','_',expName,'_',expType];

dataname    = [dataPath,'experiments','_',expName];

%% Data

di              = load(dataname);
[nprob,nsolv]   = size(di.outIts);
exs             = zeros(nprob,nsolv);

if strcmp(expType,'TIME') == true
    
    expOut = di.outTimes;
    
else
    
    expOut = di.outIts;
    
end

% Find exit criteria
for i = 1:nprob
        
    for j = 1:nsolv
        
        ng  = di.outNgs(i,j);
        
        if isnan(ng) == 1
            
            exs(i,j)    = -1;
            expOut(i,j) = -1;
            
            continue;
            
        end
        
        if ng <= convTol
        
            exs(i,j) = 1;
            
        else
            
            exs(i,j) = -1;
            
        end
        
    end    
    
end
    
if isempty(indAlg) == true
    
    indAlg = 1:nsolv;    
    
else
    
    types.colors    = types.colors(indAlg);
    types.lines     = types.lines(indAlg);
    types.markers   = types.markers(indAlg);
    
end
if isempty(legLocation) == true
    legLocation='SouthEast';
end
%perf_fnc(exs(:,indAlg),its(:,indAlg),leg(indAlg),1,types,'SouthEast',taumax);
perf_fnc(exs(:,indAlg),expOut(:,indAlg),leg(indAlg),1,types,legLocation);

%print(gcf, '-dpsc2', fullfile(figpath,'its_EX4.eps'));

fig                     = gcf;
fig.PaperPositionMode   = 'auto';
fig_pos                 = fig.PaperPosition;
fig.PaperSize           = [fig_pos(3) fig_pos(4)];

if taumax ~= 0
    
    figname = [figname,num2str(taumax),'tau.pdf'];
    
else
    
    figname = [figname,'.pdf'];
    
end

print(fig,'-dpdf',fullfile(figPath,figname));

out = fig;

end
