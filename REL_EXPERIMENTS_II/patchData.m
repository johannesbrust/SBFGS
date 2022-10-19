function patchData(dataFrom,dataTo)
% patchData: patches last column of DataTo by DataFrom
%
% INPUTS:
% dataFrom: File to load data from
% dataTo: File to load data to
%
%--------------------------------------------------------------------------
% 21/11/19, J.B., initial version

df = load(dataFrom);
dt = load(dataTo);

[np,ns] = size(df.outIts);

% Empty containers

outData = cell(np,ns);
outIts  = zeros(np,ns);
outObjs = zeros(np,ns);
outNgs  = zeros(np,ns);
outTimes= zeros(np,ns);

% To 
outData(:,1:(ns-1))     = dt.outData(:,:);
outIts(:,1:(ns-1))      = dt.outIts(:,:);
outObjs(:,1:(ns-1))     = dt.outObjs(:,:);
outNgs(:,1:(ns-1))      = dt.outNgs(:,:);
outTimes(:,1:(ns-1))    = dt.outTimes(:,:);

% From
outData(:,ns)           = df.outData(:,ns);
outIts(:,ns)            = df.outIts(:,ns);
outObjs(:,ns)           = df.outObjs(:,ns);
outNgs(:,ns)            = df.outNgs(:,ns);
outTimes(:,ns)          = df.outTimes(:,ns);

patchName = [dataTo,'_patch'];
save(patchName,'outData','outIts','outObjs','outNgs','outTimes');

end