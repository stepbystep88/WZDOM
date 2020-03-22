function y=multMatr1(A,x,flag)
%y = Ax or A'x
%
% Call: y=multMatr1(A,x,flag)
% 
%flag=0: y=A*x
%flag=1: y=A'*x

if flag,
	y=(x'*A)';
else
	y=A*x;
end