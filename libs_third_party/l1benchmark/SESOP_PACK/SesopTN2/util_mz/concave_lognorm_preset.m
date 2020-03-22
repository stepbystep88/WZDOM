function [p, h, b] = concave_lognorm_preset(b,c);
% Prepare parameters for the approximation of concave "norm" using logarithms
%
% b - second derivative at the origin
% c  - argument scaling
%
% h=  (log(1+c) - log(1+p*c)/p)
%
% We will use the calculated parameters to build the function concave_lognorm:
% phi(t)= (1/h)*(log(1+ct) -  log(1+pct)/p)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%

global Global_conc_lognorm

hinv=log(1+c);
%fprintf('p with iterations:\n')
for i=1:1000,                  % Iterative  refinement of  p 
    p=hinv * b/(c^2) +1;
    %if mod(i,50)==0, fprintf('%g ',p);end
    hinv=  (log(1+c) - log(1+p*c)/p); 
end

p=max(p,1.00001);
h = 1./ (log(1+c) - log(1+p*c)./p); 
b=h.*c^2.*(p-1);

Global_conc_lognorm.p =p;
Global_conc_lognorm.h =h;
Global_conc_lognorm.b =b;