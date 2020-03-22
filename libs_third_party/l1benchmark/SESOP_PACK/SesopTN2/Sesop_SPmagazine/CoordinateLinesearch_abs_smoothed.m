function x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps);
%Minimize in x_s        0.5w*(x_s - x0)^2 +g*x_s+ lambda0*phi(x_s,eps)
%
%      phi(x,eps) computed by abs_smoothed_eps(x,eps):
%
%           In convex case:      phi(x)= abs(x) - eps* log(abs(x)/eps +1);
%
%           In "concave" case  (when flag Global_conc_lognorm.on=1):
%
%                      phi(x)= (1/h)*(log(1+cx) -  (1/p1)*log(1+p1*cx)); 
%
%                       where     c=1/eps;
%                                        p1=Global_conc_lognorm.p;
 %                                       h  =Global_conc_lognorm.h ;
 %        
 %                       For details see:  concave_lognorm.m 
 %                                                      concave_lognorm_preset.m
 %                      


global          Global_conc_lognorm;    

p=1./eps;  % Use p instead of eps  for compatibility with old versions

% Optimality condition in x of quadratic-linear term: w*(x-x0)+g=0 gives us
 x = x0-g./w;
 
%Now we can minimize in x_s "canonic" form:   0.5(x_s - x)^2 + lambda*phi(x_s),     where
lambda=lambda0./w;




if ~isempty(Global_conc_lognorm) && Global_conc_lognorm.on, %%%%%%%%% when concave norm approximation is used: concave_lognorm()

	if Global_conc_lognorm.simple, % Simple version (nonthmooth at origin): phi(t)= (1/p)*log(1+p|t|)
		t=abs(x);
		d=p.*t-1;
      discrim=d.^2-4*p.*(lambda-t);
      tmp=double(discrim>0);  % if discriminant negative, put x_s=0;
      discrim=discrim.*tmp; 
		x_s=tmp.*sign(x).*max(0,(d+sqrt(discrim))./(2*p));

	else   % Smooth, but more complicated: phi(t)= (1/h)*(log(1+c|t|) -  (1/p)*log(1+pc|t|))

		% p,c  - argument scalings
		% h=  log(1+c) - log(1+p*c)/p)

		c=p;
		p=Global_conc_lognorm.p;
		h=Global_conc_lognorm.h ;

		r=1./lambda;      % Parameters of cubic parabola
		s= -x./lambda;
		mabs_s= -abs(s);   % s -slope of the linear term of the parabola; we need abs to treat negative slope/roots
		n=length(x);
		roots=zeros(n,3);

		for i=1:n,
			roots(i,:)=cubic_roots( r(i)*p*c^2,   c*( r(i) + p*(r(i)+c*mabs_s(i))),  r(i)+mabs_s(i)*(c+p*c) +h*c^2*(p-1),   mabs_s(i)  );
			roots(i,:) = -sign(s(i)).*roots(i,:) ;  % treat also the case when all the roots should be negative;
		end
		rroots=real(roots);    %  + 1e10*imag(roots); %this will provide min for a real root

		phi= concave_lognorm(rroots, c);
		f=repmat(r/2,1,3).*rroots.^2+repmat(s,1,3).*rroots   +  phi;
		% f=  repmat(f,1,3)
		[tmp,ind]=min(f,[],2);
		x_s=rroots(sub2ind([n,3], [1:n]',ind));
	end

else  %%%%%%%%%%%%%%%%%%%%%% When convex norm approximation is used:

	if p>=1e20,  % when we use non-smoothed abs
		
		x_s= sign(x).* max(abs(x) - lambda,0);

	else        % for smoothed abs
	
		x_s= sign(x).* soft_th_smoothed1(abs(x),lambda,p); %% Mzib
	
	end
end

%if any(imag(x_s)),keyboard;end

end




function x_s=soft_th_smoothed1(x,lambda,p)

temp=p*(x-lambda);
x_s=(temp-1+sqrt((1-temp).^2+4*p*x))/(2*p);

end




function dummy

%%%%%%%%%%%%%%%%%% Numerical testing (evaluate in the command window)
global Global_conc_lognorm
%Global_conc_lognorm.on=0;eps=1e-8;
Global_conc_lognorm.on=1;
Global_conc_lognorm.simple=1;eps=1e8;

c=20;    % argument scaling
b=100;     % second derivative at the origin (may be increased by concave_lognorm_preset to match c)
dt=1e-4;
t=[-10:dt:10];
[p1, h, b] = concave_lognorm_preset(b,c);
%[phi, dphi,d2phi] = concave_lognorm(t, c);
[phi, dphi,d2phi] = abs_smoothed_eps(t,eps);
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
eps=0.1; w=0.3;x0=-5; g=0.11;lambda0=1.7;
%eps=1e-3;lambda0=1e-3; w=0.3;x0= 4.0318e-004; g= -3.3352e-004;
f= 0.5*w*(t - x0).^2 +g*t + lambda0*abs_smoothed_eps(t,eps);
figure;plot(t,f);grid
x_s=CoordinateLinesearch_abs_smoothed(x0,g,w,lambda0,eps)

end
