n=32;
s1=randn(n,n);
h = daubcqf(4,'min');
L=2;
alpha1=mrdwt_TI2D(s1,h,L);
alpha2=randn(size(alpha1));

%c=sum(sum(alpha1.*alpha2))
c=alpha1(:)' * alpha2(:);

s2=mirdwt_TI2D(alpha2,h,L);
d=s1(:)' * s2(:);
c/d