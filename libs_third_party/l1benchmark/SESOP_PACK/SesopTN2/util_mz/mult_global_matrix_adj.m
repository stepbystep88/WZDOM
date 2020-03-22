function y = mult_global_matrix_adj(x,par)
% Compute  y=GlobalMatrix' * x
%
%Call:   y = mult_global_matrix_adj(x,par)
%
%Input:
%  x 
%  par - dummy argument
%
%Global variables:   GlobalMatrix

% Michael Zibulevsky 10.07.2008


global GlobalMatrix
y=(x'*GlobalMatrix)';
