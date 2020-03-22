function sigma_max = max_singular_value(multA,multAt,x,niter,par)
%Compute maximal singular value of matrix A, given by multiplication
%operators multA(x,par),multAt(x,par); 
%input x should be a random vector

for i=1:niter,
	x=x/norm(x);
   y=multA(x,par);
	x=multAt(y,par);
end

sigma_max = sqrt(norm(x));

return

%%%% Test
m=2;n=3;niter=50;
A=rand(m,n);
x=ones(n);
par=[];
multA=@(x,par) A*x;
multAt=@(x,par) A'*x;
svd(A)'
max_singular_value(multA,multAt,x,niter,par)