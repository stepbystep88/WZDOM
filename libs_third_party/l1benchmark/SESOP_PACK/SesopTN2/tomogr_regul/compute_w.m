% In order to approximately invert Radon transfotm, we compute
% reconstruction functional "w", which has minimal L2-norm given bluring 
% window "wind". Penalty "lam" compromises between ||w|| and the
% influence of pixels located out of the window.

%---------------------- By Michael Zibulevsky, 28.03.2006; version 14.04.2006 ------------------------------%

%profile on;compute_w;profile viewer;

%close all
clear all


par.max_iter=100; % Max iterations of Conjugate Gradients

res = 200;              % Image resolution in pixels (square image)

Nangles=400;

angles = linspace(-45,135,Nangles);   % Our fast Radon works only for angles from -45 to 135 degrees
%angles = linspace(-45,45,Nangles);;

par.k_small_bins =4; % numper of small bins per interval between pixels (used in fastradon.m)
par.flag_fastradon=1;


lam=1e-1;            %penalty par. for norm w 
mu=2;                %penalty par. for the constraint 1'A_1*w=1

k1=15;    %horizontal window size;
k2=k1;   %vertical window size;

i1=ceil(res/4);    %window corner
i2=ceil(res/4);
ii=[i1:i1+k1-1; i2:i2+k2-1]; % window indeces


wind=zeros(res);
wind(ii(1,:),ii(2,:))=1;


par.wind=wind;
par.lam=lam;
par.mu=mu;
par.angles=angles;
par.res=res;

b=myradon(wind,par);
[Nbins,Nang]=size(b);

par.b=b;
par.Nbins=Nbins;
par.Nang=Nang;

figure;
w0=zeros(size(b));
par.report=@report;                   % function for plotting CG progress results
w = CGmz(@fw,w0,mu*b,par);  % Gonjugate Gradient solution of Hw=mu*b;  fw.m computes Hw

Atw=myadjradon(w,par);

w=w/(wind(:)'*Atw(:));  % normalization 
Atw=Atw/(wind(:)'*Atw(:));

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

