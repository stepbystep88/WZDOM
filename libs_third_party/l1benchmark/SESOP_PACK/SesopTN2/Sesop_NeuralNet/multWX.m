function  vcub=multWX(vcwb,par,X)
%vcub=multA(vcwb,par):  U=W*X; vcub=[v;c;U(:);b];

% Michael Zibulevsky, 04.08.2008

M=par.nneurons;
v=vcwb(1:M);
c=vcwb(M+1);
w=vcwb(M+2 : end-M);
b=vcwb(end-M+1:end);
[N, K]=size(X);            % matrix X:  N inputs (including ones) x K examples
W=reshape(w,par.nneurons,N);
U=W*X;
vcub=[v;c;U(:);b];
