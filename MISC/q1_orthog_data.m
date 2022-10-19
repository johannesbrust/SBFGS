function [ q1, q1b, d ] = q1_orthog_data( n,r, laml, lamu, phi )
%q1_orthog_data Orthogonal rank-k perturbation to In identity matrix
%   
% Computed matrix:
%
% q1 = phi.*In + \bar{Q1} (D - Ir) \bar{Q1}',
%
% where \bar{Q1} holds orthonormal columns, \bar{Q1} \in n x r.
%
% Inputs:
% n     : Dimension
% r     : Perturbation rank
% laml  : Lower eigenvalue threshold
% lamu  : Upper eigenvalue threshold
%
% Outputs:
% q1    : Perturbed identiy
% q1b   : Eigenvectors
% d     : Eigenvalues (d+phi)
%-------------------------------------------------------------------------%
% J.B., 10/21/19
% J.B., 11/15/19, Inclusion of phi parameter

% Perturbation
d   = laml + (lamu-laml)*rand(r,1);

q1b = orth(randn(n,r));

% Full matrix
q1  = phi.*eye(n) + (q1b*diag(d))*q1b';


end

