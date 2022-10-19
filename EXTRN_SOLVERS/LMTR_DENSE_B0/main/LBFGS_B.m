function [ x,f,info ] = LBFGS_B( fnc, x0, pars )
%LBFGS_B Wrapper for the matlab interface of lbfgsb.m.

x           = x0;
[x,f,info]  = lbfgsb( fnc, pars.l, pars.u, pars) ;



end

