function outinfo=LMTR_out(numf,numg,tcpu,ex,it,tract,trrej,params,numrst)
%   LMTR_out returns output information from LMTR:
%           ex      - exitflag:
%                            1, if  norm of gradient is too small
%                           -1, if  TR radius is too small, trrad<1e-15
%                           -2, if  the max number of iterations is exceeded
%           numit   - number of succesful TR iterations
%           numf    - number of function evaluations
%           numg    - number of gradient evaluations
%           tcpu    - CPU time of algorithm execution
%           tract   - number of iterations when TR is active
%           trrej   - number of iterations when initial trial step is rejected
%           params  - input paramaters
%           numrst  - number of restarts, optional
outinfo=struct;
outinfo.numf=numf;
outinfo.numg=numg;
outinfo.tcpu=tcpu;
outinfo.ex=ex;
outinfo.numit=it;
outinfo.tract=tract;
outinfo.trrej=trrej;
outinfo.params=params;
if nargin==9
    outinfo.numrst=numrst;
end