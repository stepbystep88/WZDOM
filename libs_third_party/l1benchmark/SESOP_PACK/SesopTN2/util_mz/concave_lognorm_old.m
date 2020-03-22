function [phi, dphi,d2phi] = concave_lognorm(t, c);
% Approximation of concave "norm"  phi(t)= (1/h)*(log(1+ct) -  (1/p)*log(1+pct))
% p,c  - argument scalings
% h=  log(1+c) - log(1+p*c)/p)
%
% We pre-calculate the  parameters p,h  by   the function  concave_lognorm_preset.m
%
% Michael Zibulevsky
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Global_conc_lognorm
p=Global_conc_lognorm.p;
h=Global_conc_lognorm.h ;



u=c*abs(t);  v=p*u;
u1=u+1; v1=v+1;
phi= h*(log(u1) - (1/p)* log(v1));

if nargout>1,
    tmp1=1./u1; tmp2=1./v1;
    dphi= sign(t).*(h*c).*(tmp1-tmp2);
end

if nargout >2,
    d2phi=h*c*c*(p*tmp2.^2 - tmp1.^2);
end

return


%%%%%%%%%%%%%%%%%% Numerical testing (evaluate in the command window)
c=20;    % argument scaling
b=100;     % second derivative at the origin (may be increased by concave_lognorm_preset to match c)
dt=1e-4;
t=[-10:dt:10];
[p, h, b] = concave_lognorm_preset(b,c);
[phi, dphi,d2phi] = concave_lognorm(t, c);
figure;subplot(311);plot(t,phi);title(sprintf('phi:   b=%g  c=%g  p=%g', b, c,p));grid
dphi_num=[diff(phi)/dt  0];
subplot(312);plot(t, dphi, t,  dphi_num);title('[dphi dphinum ]');grid
d2phi_num=[diff(dphi)/dt  0];
subplot(313);plot(t, d2phi, t,  d2phi_num);title('[d2phi d2phinum ]'); grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

