function z = dwt_synthesis(v, h, levels) %mirdwt_TI2D.m
%
% 
% Wrapper for the mirdwt function of the
% Rice Wavelet Woolbox, to be used in 
% combination with mirdwt_TI.
% 
%
% These two functions apply rescalings such 
% that mirdwt_TI and midwt_TI correspond to
% multiplying by a matrix and its transpose.
%
% Written by Mario Figueiredo, 12/05/2005
%
scalefactor = 2;
[n1 n2] = size(v);
n = min(n1,n2);
t1 = v(:,1:n)*scalefactor^((levels-1));
for ll = 1:levels-1
    t2(:,(ll-1)*n*3+1:ll*n*3) = ...
    v(:,n+(ll-1)*n*3+1:n+ll*n*3)*scalefactor^(ll);
end
z = mirdwt(t1,t2,h,levels-1);

