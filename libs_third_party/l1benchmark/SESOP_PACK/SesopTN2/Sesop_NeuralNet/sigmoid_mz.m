function [f,df,d2f]=sigmoid_mz(t,eps,flagXnew)
% Sigmoid function:   f(t) = t/(1+|t|);  f_eps(t)=f(t/eps);
%Becomes more sharp  when eps --> 0
%
% Our sigmoid is a derivative of  the smooth abs. value appoximation  F= |t|-log(1+|t|);

% Michael Zibulevsky 10.10.2006;  21.07.2008; 04.08.2008
%
% Copyright (c) 2008. All rights reserved. Free for academic use. No warranty 

persistent f_old df_old d2f_old  tt tmp tmp1  flagf  flagdf  flagd2f



if nargin<3 || flagXnew
	flagf=1;flagdf=1;flagd2f=1;
end
	
p=(1./eps);

if flagf
	tt=p.*t;
	u = abs(tt);    % Our sigmoid is a derivative of  F= eps.*(u-log(1+u));
	%sign_tt=sign(tt);
	%u = sign_tt.*tt;   %abs(tt);    % Our sigmoid is a derivative of  F= eps.*(u-log(1+u));
	tmp=1./(1+u);
	%tmp=1./(1+sign_tt.*tt);
	f=tt.*tmp;
	
	f_old=f;
	%tt_old=tt;
	flagf=0;
else
	f=f_old;
	%tt=tt_old;
end


if nargout > 1,      % First derivative
	if flagdf
		tmp1=tmp.*tmp;
		df= p.*tmp1;
		
		df_old=df;
		flagdf=0;
	else
        df=df_old;
	end		
end

if nargout > 2,      % Second derivative
		if flagd2f
			%d2f= -p.*(tt./u).*tmp1.*tmp;
			%d2f= (-2*p.^2).*sign_tt.*tmp1.*tmp;
			d2f= (-2*p.^2).*sign(tt).*tmp1.*tmp;
			d2f_old=d2f;
			flagd2f=0;
		else
			d2f=d2f_old;
		end
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             TEST:  execute this part as a script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a=5;
da=0.001;
x=[-a:da:a]';
n=length(x);

[phi,dphi,d2phi]= sigmoid_mz(x, 0.7);
figure;plot(x,[phi,dphi,d2phi]);grid
figure;plot(x(1:n-1),[(phi(2:n)-phi(1:n-1))/da dphi(1:n-1)]);grid
figure;plot(x(1:n-1),[(dphi(2:n)-dphi(1:n-1))/da d2phi(1:n-1)]);grid

% figure;plot(x(1:n-1),[d2phi(1:n-1)]);grid
% figure;plot(x(1:n-1),[dphi(1:n-1)]);grid
% figure;plot(x(1:n-1),(phi(2:n)-phi(1:n-1))/da);grid
% figure;plot(x(1:n-1),(dphi(2:n)-dphi(1:n-1))/da);grid
