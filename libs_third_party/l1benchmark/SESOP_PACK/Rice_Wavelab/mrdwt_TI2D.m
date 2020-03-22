function z = mrdwt_TI(v, h, levels)
%
% 
% Wrapper for the mrdwt function of the
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
[m n] = size(v);
[temp1,temp2,lev] = mrdwt(v,h,levels-1);
temp1 = temp1*scalefactor^(-(levels-1));
for ll = 1:levels-1
     temp2(:,(ll-1)*n*3+1:ll*n*3) = ...
     temp2(:,(ll-1)*n*3+1:ll*n*3)*scalefactor^(-ll);
end
z = [temp1 temp2];

