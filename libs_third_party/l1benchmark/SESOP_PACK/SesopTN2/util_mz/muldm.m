function B=muldm(d,A) 
% Compute efficiently  B= diag(d)*A
%
% Call:
%B=muldm(d,A)
%

d=d(:);
[M,N]=size(A);
n=length(d);
if n ~= M, error('Wrong matrix sizes');end

B=A;
for j=1:N
B(:,j)=d.*B(:,j);
end