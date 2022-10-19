function [ H ] = RosenbrockHessian( x, params )
%ROSENBROCKHESSIAN Computes the 2nd derivative matrix of the Rosenbrock
%function.


a= params{1};
b= params{2};

twoM    = size(x,1);
M       = twoM/2;

H       = zeros(twoM,twoM);

for i = 1:M
    
   H(2*i,2*i)       = 2*b;
   H(2*i,2*i-1)     = -4*b*x(2*i-1);
   H(2*i-1,2*i)     = -4*b*x(2*i-1);
   H(2*i-1,2*i-1)   = -4*b*( x(2*i) - x(2*i-1)^2 + x(2*i-1)*( x(2*i) - 2*x(2*i-1) ) ) + 2;

end

