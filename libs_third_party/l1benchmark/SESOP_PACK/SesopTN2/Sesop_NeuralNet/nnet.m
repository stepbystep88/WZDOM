function y=nnet(vcwb,par,X)
%Simgle layer Neural Net:  y=v'*sigmoid(W*X+b*1')+c
%
% Call:  y = nnet(vcwb,par)   
% Fields in use: par.nneurons,  par.eps_sigmoid

% Michael Zibulevsky 24.07.2008; 04.08.2008

[N, K]=size(X);    % matrix X:  N inputs (including ones) x K examples
M=par.nneurons;

v=vcwb(1:M);
c=vcwb(M+1);
w=vcwb(M+2 : end-M);
b=vcwb(end-M+1:end);

W=reshape(w,M,N);
Phi=sigmoid_mz(W*X+b*ones(1,K), par.eps_sigmoid);
y=v' * Phi +c ;  % output of nnet
