%------------------------ test_LIBSVM_fun --------------------------------%
%
% test_LIBSVM_fun.m is a test for the LIBSVM function and gradient
% computation.
%-------------------------------------------------------------------------%
% 10/04/19, J.B.
% 10/15/19, J.B., testing the vectorized version

global LAM;

LAM = 1e-3;

n       = 1e2;
m       = 1e1;
  
X       = randn(m,n);

ry      = rand(m,1);

y       = ones(m,1);

% Definition of labels
y(ry < 0.5)  = -y(ry < 0.5);

w       = randn(n,1);

f       = LIBSVM_fun(w,y,X);

[f1,g]  = LIBSVM_fun(w,y,X);

e1      = f-f1;

% Wrapper function test

fw      = @(z)( LIBSVM_fun(z,y,X) );

ffw         = fw(w);
[ffw1,gw]   = fw(w);

e2 = f-ffw1;
  
% Vectorized calls

fv          = LIBSVM_fun_VEC(w,y,X);

[f1v,gv]    = LIBSVM_fun_VEC(w,y,X);

e3          = f-f1v;
e4          = norm(g-gv);