% This is a script to minimize Rozenbrock function using the limited
% memory trust-region algorithm EIG(inf,2). 
%
% For details, see:
% O.Burdakov, L. Gong, Y. Yuan, S. Zikrin, On Efficiently Combining Limited
% Memory and Trust-Region Techniques, technical report LiTH-MAT-R-2013/13-SE,
% Department of Mathematics, Linköping University, 2013.
% http://liu.diva-portal.org/smash/record.jsf?pid=diva2%3A667359

parentpath = [cd(cd('..')) '/main'];
addpath(parentpath)

% Set input parameters:
params=struct;
params.m = 5;  % number of L-BFGS updates
params.gtol = 1e-5;  % exit if ||g||_2<gtol*max(1,||x||_2)
params.ranktol = 1e-7;  % tolerance for establishing rank of V
params.dflag = 0;  % display parameter, 1 if to display information per iteration
params.trtol = 0.1;  % exit MS algorithm if abs(||s||-delta)<trtol*delta
params.ftol=1e-11;  % tolerance on relative function reduction

% Test on Rosenbrock function
x0=-ones(1000,1);
%[x,f,outinfo] = LMTR_EIG_inf_2(@rosen,x0,params);
[x,f,outinfo] = LMTR_EIG_MS_2_2(@rosen,x0,params);
