function A=spdiag(d)
%Diagonal matrix in sparse format
n=length(d);
ind=[1:n];
A=sparse(ind,ind,d,n,n);
