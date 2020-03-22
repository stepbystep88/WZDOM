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
        % phi(t)= (1/h)*(log(1+ct) -  (1/p)*log(1+pct))
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

    return
    
end



%%%%%%%%%%%%%%%%%%%%%% When convex norm approximation is used:

%mask=double(x>=0);
% x_s=soft_th_smoothed1(x.*mask,lambda,p)-soft_th_smoothed1(-x.*(1-mask),lambda,p);

x_s= sign(x).* soft_th_smoothed1(abs(x),lambda,p); %% Mzib






function x_s=soft_th_smoothed1(x,lambda,p)

temp=p*(x-lambda);
x_s=(temp-1+sqrt((1-temp).^2+4*p*x))/(2*p);

