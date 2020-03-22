function  vcwb=multUXt(vcu,par,X)
%Adjoint operator: W=U*X'; vcwb=[v;c;W(:);b];

M=par.nneurons;
K=par.Ktrain_samples;

indv=[1:M];
indc=M+1;
%indU=[M+2: M+1+M*K];
indb=[ (M*(K+1) +2) : (M*(K+2) +1)];

v=vcu(indv);
c=vcu(indc);
b=vcu(indb);
u=vcu([(M+2): (M+1+M*K)]);
U=reshape(u, M, K);

W=U*X';
vcwb=[v;c;W(:);b];

