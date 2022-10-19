function [ex,numf,numg,numit,tcpu,tract,numrst]=runAlgorithm(alg,x0,params,numRuns)
tcpu=zeros(numRuns,1);
for i=1:numRuns
    [x,f,outinfo]=alg(@cutest_fun,x0,params);
    %[~,~,outinfo]=alg(@cuter_fun,x0,params);
    tcpu(i)=outinfo.tcpu;
end
numf=outinfo.numf;
numg=outinfo.numg;
ex=outinfo.ex;
numit=outinfo.numit;   
tract=outinfo.tract;
if nargout==7
    numrst=outinfo.numrst;
end