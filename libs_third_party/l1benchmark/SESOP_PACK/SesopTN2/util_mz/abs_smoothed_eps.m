function [f,g,h]=abs_smoothed_eps(x,eps) 
%Convex:   |x| - eps* log(|x|/eps +1);  "Concave":  (1/h)*(log(1+c|x|) -  (1/p1)*log(1+p1*c|x|)); 
%
%           In "concave" case  (when flag Global_conc_lognorm.on=1):
%
%                      phi(x)= (1/h)*(log(1+c|x|) -  (1/p1)*log(1+p1*c|x|)); 
%
%                       where     c=1/eps;
%                                        p1=Global_conc_lognorm.p;
 %                                       h  =Global_conc_lognorm.h ;
 %        
 %                       For details see:  concave_lognorm.m 
 %                                                      concave_lognorm_preset.m
 %                      
 %  Global variable:       Global_conc_lognorm (used only in "concave" case) 
 %
 % Michael Zibulevsky 16.07.2008 
 
% MZ:  uncomment for  the  penalty for negativity 


%Quadratic function, uncomment for debugging:
% f=0.5*x.^2;  g=x;  h=ones(size(x));  return

p=1./eps;  % Use p instead of eps  for compatibility with old versions


%%%%%%%%%%%%%%%%%%%%%%%% %%   Mzib: for Lp norm approximation, p<1
global          Global_conc_lognorm;                   
if  ~isempty(Global_conc_lognorm) && Global_conc_lognorm.on,
    
    if nargout==1,              f = concave_lognorm(x, p);  % p is used as c inside  concave_lognorm()
    elseif nargout==2, [f,g] = concave_lognorm(x, p);
    else                              [f,g,h] = concave_lognorm(x, p);
    end

    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if p>=1e20,  % Non-smoothed abs

	f=abs(x);

	if nargout>1
		g=sign(x);
	end

	if nargout>2
		h=0*x;
	end

else  %  smoothed abs

	temp=p*abs(x);
	f=(1/p)*(temp-log(temp+1));
	%f=(temp-log(temp+1))/p - x; % MZ: penalty for negativity

	if nargout>1
		g=temp.*sign(x)./(1+temp);
		%g=temp.*sign(x)./(1+temp) - 1; % MZ: penalty for negativity
	end
	if nargout>2
		h=p./( (1+temp).^2 );
	end

end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TEST   convex case:
global Global_conc_lognorm;
Global_conc_lognorm.on  =  0;
x=[-10:0.1:10]'; p=10; y=abs_smoothed(x,p);figure;plot(x,y);grid

% TEST   concave case:
global Global_conc_lognorm;
Global_conc_lognorm.on  =  1;
b    =   10; c=20; 
[p_tmp, h_tmp, b_tmp] = concave_lognorm_preset(b,  c)
x=[-10:0.001:10]'; y=abs_smoothed(x,c);figure;plot(x,y);grid

 %f_1 = abs(x) - log(abs(x) +1);      f_eps = eps* f_1(x/eps)