function [phi, dphi,d2phi] = concave_lognorm(t, c);
% Approximation of concave "norm"  phi(t)= (1/h)*(log(1+c|t|) -  (1/p)*log(1+pc|t|))
% p,c  - argument scalings
% h=  log(1+c) - log(1+p*c)/p)
%
% We pre-calculate the  parameters p,h  by   the function  concave_lognorm_preset.m
%
% Michael Zibulevsky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Global_conc_lognorm

u=c*abs(t); u1=u+1;

if Global_conc_lognorm.simple, % Simple version (nonthmooth at origin): phi(t)= (1/c)*log(1+c|t|) 
	phi=(1/c)*log(u1);
	if nargout>1,
		tmp1=1./u1; 
		dphi= sign(t).*tmp1;
	end

	if nargout >2,
		d2phi= (-c)*tmp1.^2;
	end

else          % Smooth, but more complicated: phi(t)= (1/h)*(log(1+c|t|) -  (1/p)*log(1+pc|t|))

	p=Global_conc_lognorm.p;
	h=Global_conc_lognorm.h ;

	v=p*u; v1=v+1;
	phi= h*(log(u1) - (1/p)* log(v1));

	if nargout>1,
		tmp1=1./u1; tmp2=1./v1;
		dphi= sign(t).*(h*c).*(tmp1-tmp2);
	end

	if nargout >2,
		d2phi=(h*c*c)*(p*tmp2.^2 - tmp1.^2);
	end

end

return


%%%%%%%%%%%%%%%%%% Numerical testing (evaluate in the command window)
global Global_conc_lognorm
Global_conc_lognorm.on=1;
Global_conc_lognorm.simple=1;

c=20;    % argument scaling
b=100;     % second derivative at the origin (may be increased by concave_lognorm_preset to match c)
dt=1e-4;
t=[-10:dt:10];
[p1, h, b] = concave_lognorm_preset(b,c);
[phi, dphi,d2phi] = concave_lognorm(t, c);
figure;subplot(311);plot(t,phi);title(sprintf('phi:   b=%g  c=%g  p=%g', b, c,p1));grid
dphi_num=[diff(phi)/dt  0];
subplot(312);plot(t, dphi, t,  dphi_num);title('[dphi dphinum ]');grid
d2phi_num=[diff(dphi)/dt  0];
subplot(313);plot(t, d2phi, t,  d2phi_num);title('[d2phi d2phinum ]'); grid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test coordinate minimization 
% x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps);
%
% Minimize in x_s        0.5w*(x_s - x0)^2 +g*x_s+ lambda0*phi(x_s,eps)
%
eps=0.1;
w=0.3;x0=-5; g=0.11;lambda0=1.7;
f= 0.5*w*(t - x0).^2 +g*t + lambda0*abs_smoothed_eps(t,eps);
figure;plot(t,f);grid
x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps)






