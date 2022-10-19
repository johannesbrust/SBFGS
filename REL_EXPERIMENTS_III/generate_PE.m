function [ yh, g, A] = generate_PE( n, xl, xu, yl, yu, a, b, sig )
%generate_PE: Poisson Data for PDE constrained optimization
%--------------------------------------------------------------------------
%
% Function to generate data for PDE constrained optimization. Based on
% Data required for solve in a Poisson problem:
%
% Laplace(yh)   = f(x,y)    (Interior)
%           yh  = g         (Boundary)
%
% For this problem: f(x,y) = exp(-1/sig((x-a)^2 + (y-b)^2)), with
% box domain.
%
% INPUTS:
% n:=   Mesh size
% xl:=  Lower x domain
% xu:=  Upper x domain
% yl:=  Lower y domain
% yu:=  Upper y domain
% a:=   Center for x direction
% b:=   Center for y direction
% sig:= Standard deviation
% OUTPUTS:
% yh:= Computed finite-difference solution
% g:= Boundary values. Matrix g = [L T R B]/[Left,Top,Right,Bottom]
% A:= Finite-difference matrix
%--------------------------------------------------------------------------
% 11/04/19, J.B.

    nm2   = (n-2);
    nm2sq = nm2*nm2;

    % 2D Mesh
    x     = linspace(xl,xu,n);
    y     = linspace(yl,yu,n);
    h     = 1/(n*n);
    [X,Y] = meshgrid(x,y);

    % Right-hand sides
    expXYs      = exp(-((X-a).^2-(Y-b).^2)./sig);
    fRhs        = ((4.*((X-a).^2+(Y-b).^2)./sig-4)./sig).*expXYs;

    GRhs        = -h.*fRhs((2:end-1),(2:end-1));

    GRhs(1,:)   = GRhs(1,:) + expXYs(1,(2:end-1));
    GRhs(:,1)   = GRhs(:,1) + expXYs((2:end-1),1);
    GRhs(end,:) = GRhs(end,:) + expXYs(end,(2:end-1));
    GRhs(:,end) = GRhs(:,end) + expXYs((2:end-1),end);

    gRhs        = reshape(GRhs,nm2sq,1);

    % Poisson matrix
    en2     = ones(nm2sq,1);
    en      = ones(nm2,1);
    Ap      = spdiags([-en 4*en -en],-1:1,nm2,nm2);
    ApL     = kron(speye(nm2),Ap);

    A       = ApL + spdiags(-en2,nm2,nm2sq,nm2sq) + spdiags(-en2,-nm2,nm2sq,nm2sq);

    yh      = A\gRhs;

    % Boundary data
    g       = zeros(nm2,4);

    g(:,1)  = expXYs((2:end-1),1);
    g(:,2)  = expXYs(1,(2:end-1));
    g(:,3)  = expXYs((2:end-1),end);
    g(:,4)  = expXYs(end,(2:end-1));

end

