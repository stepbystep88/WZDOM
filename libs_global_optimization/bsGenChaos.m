%  10 Chaotic maps for meta-heuristics                                               %
%                                                                                    %
%  Developed in MATLAB R2014a                                                        %
%                                                                                    %
%  Author and programmer: Seyedali Mirjalili                                         %
%                                                                                    %
%         e-Mail: ali.mirjalili@gmail.com                                            %
%                 seyedali.mirjalili@griffithuni.edu.au                              %
%                                                                                    %
%       Homepage: http://www.alimirjalili.com                                        %
%                                                                                    %                                          %
% You can simply use any of these chaotic maps when you need a random                %
% number in [0 1] in your meta-heursitic. Details can be found in the                %
% following reference:                                                               %
% 
% S. Saremi, S. Mirjalili, A. Lewis, Biogeography-based optimisation with chaos, 
% Neural Computing and Applications, In press, 2014, Springer,
% DOI: http://dx.doi.org/10.1007/s00521-014-1597-x
function O=bsGenChaos(Index,Initial_Value,Max_iter)
Max_iter = Max_iter - 1;
x(1)=Initial_Value;
switch Index
%Chebyshev map
    case 1
for i=1:Max_iter
    x(i+1)=cos(i*acos(x(i)));
    G(i)=((x(i)+1)*100)/2;
end
    case 2
%Circle map
a=0.5;
b=0.2;
for i=1:Max_iter
    x(i+1)=mod(x(i)+b-(a/(2*pi))*sin(2*pi*x(i)),1);
    G(i)=x(i)*100;
end
    case 3
%Gauss/mouse map
for i=1:Max_iter
    if x(i)==0
        x(i+1)=0;
    else
        x(i+1)=1/x(i) - floor(1/x(i));
    end
    G(i)=x(i)*100;
end
    case 4
%Iterative map
a=0.7;
for i=1:Max_iter
    x(i+1)=sin((a*pi)/x(i));
    G(i)=((x(i)+1)*100)/2;
end
%normalize it from [-1 1] to [0 1]
a=-1; b=1; c=0; d=1;
x=((x-a)*(d-c))/(b-a);
    case 5
%Logistic map
a=4;
for i=1:Max_iter
    x(i+1)=a*x(i)*(1-x(i));
    G(i)=x(i)*100;
end
    case 6
%Piecewise map
P=0.4;
for i=1:Max_iter
    if x(i)>=0 && x(i)<P
        x(i+1)=x(i)/P;
    end
    if x(i)>=P && x(i)<0.5
        x(i+1)=(x(i)-P)/(0.5-P);
    end
    if x(i)>=0.5 && x(i)<1-P
        x(i+1)=(1-P-x(i))/(0.5-P);
    end
    if x(i)>=1-P && x(i)<1
        x(i+1)=(1-x(i))/P;
    end    
    G(i)=x(i)*100;
end
    case 7
%Sine map
for i=1:Max_iter
     x(i+1) = sin(pi*x(i));
     G(i)=(x(i))*100;
end
 
    case 8
 %Singer map 
 u=1.07;
 for i=1:Max_iter
     x(i+1) = u*(7.86*x(i)-23.31*(x(i)^2)+28.75*(x(i)^3)-13.302875*(x(i)^4));
     G(i)=(x(i))*100;
 end
 
    case 9
%Sinusoidal map
 for i=1:Max_iter
     x(i+1) = 2.3*x(i)^2*sin(pi*x(i));
     G(i)=(x(i))*100;
 end
 
    case 10
 %Tent map
 x(1)=0.6;
 for i=1:Max_iter
     if x(i)<0.7
         x(i+1)=x(i)/0.7;
     end
     if x(i)>=0.7
         x(i+1)=(10/3)*(1-x(i));
     end
     G(i)=(x(i))*100;
 end
end
O=x;